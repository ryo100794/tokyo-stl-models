#!/usr/bin/env python3
"""Convert PLATEAU CityGML building posList geometry to STL."""

from __future__ import annotations

import argparse
import math
import struct
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

import requests
from lxml import etree
from pyproj import CRS, Transformer
from tqdm import tqdm

from model_coordinates import ModelCoordinates


POSLIST_TAG = "{*}posList"
WGS84 = "EPSG:4326"
TOKYO_METERS = "EPSG:6677"


@dataclass
class BBox:
    west: float
    east: float
    south: float
    north: float

    def contains(self, lon: float, lat: float) -> bool:
        return self.west <= lon <= self.east and self.south <= lat <= self.north


@dataclass
class Stats:
    scanned_files: int = 0
    used_files: int = 0
    scanned_polygons: int = 0
    used_polygons: int = 0
    vertices: int = 0
    faces: int = 0


@dataclass
class Bounds:
    min_x: float = math.inf
    min_y: float = math.inf
    min_z: float = math.inf
    max_x: float = -math.inf
    max_y: float = -math.inf
    max_z: float = -math.inf

    def update(self, point: tuple[float, float, float]) -> None:
        x, y, z = point
        self.min_x = min(self.min_x, x)
        self.min_y = min(self.min_y, y)
        self.min_z = min(self.min_z, z)
        self.max_x = max(self.max_x, x)
        self.max_y = max(self.max_y, y)
        self.max_z = max(self.max_z, z)

    def center(self) -> tuple[float, float, float]:
        return (
            (self.min_x + self.max_x) * 0.5,
            (self.min_y + self.max_y) * 0.5,
            (self.min_z + self.max_z) * 0.5,
        )


@dataclass
class Ring:
    points_m: list[tuple[float, float, float]]
    centroid_lon: float
    centroid_lat: float


@dataclass
class TransformSet:
    ll_to_m: Transformer
    m_to_ll: Transformer


def make_transformers() -> TransformSet:
    return TransformSet(
        ll_to_m=Transformer.from_crs(WGS84, TOKYO_METERS, always_xy=True),
        m_to_ll=Transformer.from_crs(TOKYO_METERS, WGS84, always_xy=True),
    )


def command_download(args: argparse.Namespace) -> int:
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(args.url, stream=True, timeout=60) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))
        with out_path.open("wb") as fh:
            with tqdm(total=total, unit="B", unit_scale=True, desc="download") as bar:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if not chunk:
                        continue
                    fh.write(chunk)
                    bar.update(len(chunk))

    print(f"downloaded: {out_path}")
    return 0


def command_unzip(args: argparse.Namespace) -> int:
    zip_path = Path(args.zip)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zf:
        members = zf.infolist()
        for member in tqdm(members, desc="unzip", unit="file"):
            zf.extract(member, out_dir)

    print(f"unzipped: {zip_path} -> {out_dir}")
    return 0


def is_citygml_file(path: Path) -> bool:
    suffix = path.suffix.lower()
    return suffix in {".gml", ".xml", ".citygml"}


def is_building_citygml(path: Path) -> bool:
    lower_parts = [part.lower() for part in path.parts]
    lower_name = path.name.lower()
    if "bldg" in lower_name:
        return True
    for idx, part in enumerate(lower_parts):
        if part == "udx" and idx + 1 < len(lower_parts) and lower_parts[idx + 1] == "bldg":
            return True
    return False


def find_building_files(root: Path, max_files: Optional[int], file_offset: int = 0) -> list[Path]:
    files: list[Path] = []
    matched = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if not is_citygml_file(path):
            continue
        if not is_building_citygml(path):
            continue
        if matched < file_offset:
            matched += 1
            continue
        matched += 1
        files.append(path)
        if max_files is not None and len(files) >= max_files:
            break
    return files


def filter_files_by_size(files: list[Path], max_size_mb: Optional[float]) -> list[Path]:
    if max_size_mb is None:
        return files
    max_size = int(max_size_mb * 1024 * 1024)
    out: list[Path] = []
    for path in files:
        try:
            if path.stat().st_size <= max_size:
                out.append(path)
            else:
                print(f"warning: skipped large file: {path}", file=sys.stderr)
        except OSError as exc:
            print(f"warning: skipped unreadable file: {path}: {exc}", file=sys.stderr)
    return out


def infer_dimension(value_count: int, attr_dim: Optional[str]) -> int:
    if attr_dim:
        try:
            dim = int(attr_dim)
            if dim in (2, 3):
                return dim
        except ValueError:
            pass
    if value_count % 3 == 0:
        return 3
    if value_count % 2 == 0:
        return 2
    return 3


def is_lat(value: float) -> bool:
    return 20.0 <= value <= 50.0


def is_lon(value: float) -> bool:
    return 120.0 <= value <= 160.0


def classify_xy(a: float, b: float) -> str:
    if is_lat(a) and is_lon(b):
        return "latlon"
    if is_lon(a) and is_lat(b):
        return "lonlat"
    return "plane"


def parse_poslist_text(text: str, attr_dim: Optional[str], transforms: TransformSet,
                       srs_name: Optional[str] = None) -> Optional[Ring]:
    parts = text.split()
    if len(parts) < 6:
        return None

    try:
        values = [float(item) for item in parts]
    except ValueError:
        return None

    dim = infer_dimension(len(values), attr_dim)
    if len(values) < dim * 3:
        return None

    usable = len(values) - (len(values) % dim)
    coords = [values[i : i + dim] for i in range(0, usable, dim)]
    if len(coords) < 3:
        return None

    first = coords[0]
    mode = classify_xy(first[0], first[1])
    if srs_name:
        # EPSG:6677's native GML axis order is northing, easting. The public
        # model API and always_xy pyproj calls instead use easting, northing.
        try:
            if CRS.from_user_input(srs_name).to_epsg() == 6677:
                mode = "north_east"
        except Exception as exc:
            raise ValueError(f"Unrecognized declared CityGML CRS: {srs_name}") from exc
    points_m: list[tuple[float, float, float]] = []
    lon_lat_items: list[tuple[float, float]] = []

    for coord in coords:
        a = coord[0]
        b = coord[1]
        z = coord[2] if dim >= 3 else 0.0

        if mode == "latlon":
            lat = a
            lon = b
            x, y = transforms.ll_to_m.transform(lon, lat)
        elif mode == "lonlat":
            lon = a
            lat = b
            x, y = transforms.ll_to_m.transform(lon, lat)
        elif mode == "north_east":
            x, y = b, a
            lon, lat = transforms.m_to_ll.transform(x, y)
        else:
            x = a
            y = b
            lon, lat = transforms.m_to_ll.transform(x, y)

        if not all(math.isfinite(v) for v in (x, y, z, lon, lat)):
            return None
        points_m.append((x, y, z))
        lon_lat_items.append((lon, lat))

    if len(points_m) >= 2 and near_same_point(points_m[0], points_m[-1]):
        points_m.pop()
        lon_lat_items.pop()

    if len(points_m) < 3:
        return None

    centroid_lon = sum(item[0] for item in lon_lat_items) / len(lon_lat_items)
    centroid_lat = sum(item[1] for item in lon_lat_items) / len(lon_lat_items)
    return Ring(points_m=points_m, centroid_lon=centroid_lon, centroid_lat=centroid_lat)


def near_same_point(a: tuple[float, float, float], b: tuple[float, float, float]) -> bool:
    return abs(a[0] - b[0]) < 1.0e-8 and abs(a[1] - b[1]) < 1.0e-8 and abs(a[2] - b[2]) < 1.0e-8


def local_name(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", 1)[1]
    return tag


def elem_matches_lod(elem, lod: str) -> bool:
    if lod == "any":
        return True
    prefix = f"lod{lod}"
    parent = elem
    while parent is not None:
        if local_name(parent.tag).lower().startswith(prefix):
            return True
        parent = parent.getparent()
    return False


def iter_rings(path: Path, transforms: TransformSet, lod: str = "any") -> Iterator[Ring]:
    context = etree.iterparse(
        str(path),
        events=("end",),
        tag=POSLIST_TAG,
        recover=True,
        huge_tree=True,
    )

    for _event, elem in context:
        if elem_matches_lod(elem, lod):
            text = elem.text or ""
            parent = elem
            srs_name, dimension = None, None
            while parent is not None:
                srs_name = srs_name or parent.get("srsName")
                dimension = dimension or parent.get("srsDimension")
                parent = parent.getparent()
            if srs_name is None:
                envelope = elem.getroottree().getroot().find("./{*}boundedBy/{*}Envelope")
                if envelope is not None:
                    srs_name = envelope.get("srsName")
                    dimension = dimension or envelope.get("srsDimension")
            ring = parse_poslist_text(text, dimension, transforms, srs_name)
            if ring is not None:
                yield ring
        elem.clear()
        parent = elem.getparent()
        while parent is not None and elem.getprevious() is not None:
            del parent[0]


def append_ring_mesh(ring: Ring, vertices: list[tuple[float, float, float]], faces: list[tuple[int, int, int]]) -> None:
    base = len(vertices)
    vertices.extend(ring.points_m)
    count = len(ring.points_m)
    for idx in range(1, count - 1):
        faces.append((base, base + idx, base + idx + 1))


def ring_triangle_count(ring: Ring) -> int:
    return max(0, len(ring.points_m) - 2)


def scale_point(
    point: tuple[float, float, float],
    center: tuple[float, float, float],
    scale: float,
) -> tuple[float, float, float]:
    return (
        (point[0] - center[0]) * scale,
        (point[1] - center[1]) * scale,
        (point[2] - center[2]) * scale,
    )


def normal_for_triangle(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> tuple[float, float, float]:
    ux = b[0] - a[0]
    uy = b[1] - a[1]
    uz = b[2] - a[2]
    vx = c[0] - a[0]
    vy = c[1] - a[1]
    vz = c[2] - a[2]
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    length = math.sqrt(nx * nx + ny * ny + nz * nz)
    if length <= 0.0:
        return (0.0, 0.0, 0.0)
    return (nx / length, ny / length, nz / length)


def write_binary_stl_triangle(
    fh,
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> None:
    n = normal_for_triangle(a, b, c)
    fh.write(struct.pack("<12fH", n[0], n[1], n[2], a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], 0))


def scaled_bounds(bounds: Bounds, center: tuple[float, float, float], scale: float) -> Bounds:
    out = Bounds()
    for x in (bounds.min_x, bounds.max_x):
        for y in (bounds.min_y, bounds.max_y):
            for z in (bounds.min_z, bounds.max_z):
                out.update(scale_point((x, y, z), center, scale))
    return out


def coordinates_for_bounds(args: argparse.Namespace, bounds: Bounds) -> ModelCoordinates:
    """Use one explicit origin across batches, or retain legacy auto-centering."""
    if getattr(args, "origin_json", None):
        coordinates = ModelCoordinates.read(args.origin_json, mm_per_m=args.scale,
            vertical_exaggeration=getattr(args, "vertical_exaggeration", 1.0),
            z_offset_mm=getattr(args, "z_offset", 0.0))
        if CRS.from_user_input(coordinates.crs).to_epsg() != 6677:
            raise ValueError("PLATEAU output origin must use EPSG:6677 projected metres")
        return coordinates
    return ModelCoordinates(*bounds.center(), mm_per_m=args.scale,
        vertical_exaggeration=getattr(args, "vertical_exaggeration", 1.0),
        z_offset_mm=getattr(args, "z_offset", 0.0))


def transformed_bounds(bounds: Bounds, coordinates: ModelCoordinates) -> Bounds:
    result = Bounds()
    for x in (bounds.min_x, bounds.max_x):
        for y in (bounds.min_y, bounds.max_y):
            for z in (bounds.min_z, bounds.max_z):
                result.update(coordinates.point((x, y, z)))
    return result


def base_plate_triangles() -> list[tuple[int, int, int]]:
    return [
        (0, 1, 2),
        (0, 2, 3),
        (4, 6, 5),
        (4, 7, 6),
        (0, 4, 5),
        (0, 5, 1),
        (1, 5, 6),
        (1, 6, 2),
        (2, 6, 7),
        (2, 7, 3),
        (3, 7, 4),
        (3, 4, 0),
    ]


def make_base_plate_vertices(bounds: Bounds, thickness: float) -> list[tuple[float, float, float]]:
    min_x = bounds.min_x
    max_x = bounds.max_x
    min_y = bounds.min_y
    max_y = bounds.max_y
    top_z = bounds.min_z
    bottom_z = bounds.min_z - thickness
    pad_x = max((max_x - min_x) * 0.02, 1.0)
    pad_y = max((max_y - min_y) * 0.02, 1.0)
    min_x -= pad_x
    max_x += pad_x
    min_y -= pad_y
    max_y += pad_y
    return [
        (min_x, min_y, top_z),
        (max_x, min_y, top_z),
        (max_x, max_y, top_z),
        (min_x, max_y, top_z),
        (min_x, min_y, bottom_z),
        (max_x, min_y, bottom_z),
        (max_x, max_y, bottom_z),
        (min_x, max_y, bottom_z),
    ]


def scan_stl_inputs(
    files: list[Path],
    bbox: BBox,
    transforms: TransformSet,
    lod: str,
) -> tuple[Stats, Bounds, int]:
    stats = Stats(scanned_files=len(files))
    bounds = Bounds()
    face_count = 0

    for path in tqdm(files, desc="scan", unit="file"):
        used_this_file = False
        try:
            for ring in iter_rings(path, transforms, lod):
                stats.scanned_polygons += 1
                if not bbox.contains(ring.centroid_lon, ring.centroid_lat):
                    continue
                tri_count = ring_triangle_count(ring)
                if tri_count <= 0:
                    continue
                for point in ring.points_m:
                    bounds.update(point)
                stats.used_polygons += 1
                stats.vertices += len(ring.points_m)
                face_count += tri_count
                used_this_file = True
        except Exception as exc:
            print(f"warning: skipped CityGML file: {path}: {exc}", file=sys.stderr)
            continue
        if used_this_file:
            stats.used_files += 1

    stats.faces = face_count
    return stats, bounds, face_count


def write_streaming_stl(
    files: list[Path],
    out_path: Path,
    bbox: BBox,
    transforms: TransformSet,
    bounds: Bounds,
    face_count: int,
    scale: float,
    base_thickness: float,
    lod: str,
    coordinates: Optional[ModelCoordinates] = None,
) -> tuple[int, int]:
    center = bounds.center()
    coordinates = coordinates or ModelCoordinates(*center, mm_per_m=scale)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    base_faces = 12 if base_thickness > 0.0 else 0
    total_faces = face_count + base_faces
    if total_faces > 0xFFFFFFFF:
        raise RuntimeError("too many STL faces for binary STL")

    with out_path.open("wb") as fh:
        header = b"PLATEAU CityGML to STL".ljust(80, b" ")
        fh.write(header)
        fh.write(struct.pack("<I", total_faces))

        for path in tqdm(files, desc="write", unit="file"):
            try:
                for ring in iter_rings(path, transforms, lod):
                    if not bbox.contains(ring.centroid_lon, ring.centroid_lat):
                        continue
                    points = [coordinates.point(point) for point in ring.points_m]
                    for idx in range(1, len(points) - 1):
                        write_binary_stl_triangle(fh, points[0], points[idx], points[idx + 1])
            except Exception as exc:
                print(f"warning: skipped CityGML file while writing: {path}: {exc}", file=sys.stderr)
                continue

        if base_thickness > 0.0:
            sb = transformed_bounds(bounds, coordinates)
            base_vertices = make_base_plate_vertices(sb, base_thickness)
            for tri in base_plate_triangles():
                write_binary_stl_triangle(fh, base_vertices[tri[0]], base_vertices[tri[1]], base_vertices[tri[2]])

    base_vertices_count = 8 if base_thickness > 0.0 else 0
    return base_vertices_count, base_faces


def scaled_mesh_data(
    vertices_m: list[tuple[float, float, float]],
    faces: list[tuple[int, int, int]],
    scale: float,
    coordinates: Optional[ModelCoordinates] = None,
) -> tuple[np.ndarray, np.ndarray]:
    import numpy as np

    vertices = np.asarray(vertices_m, dtype=np.float64)
    center = (vertices.min(axis=0) + vertices.max(axis=0)) * 0.5
    coordinates = coordinates or ModelCoordinates(*center, mm_per_m=scale)
    scaled = coordinates.points(vertices)
    return scaled, np.asarray(faces, dtype=np.int64)


def add_base_plate(
    vertices: np.ndarray,
    faces: np.ndarray,
    thickness: float,
) -> tuple[np.ndarray, np.ndarray]:
    import numpy as np

    if thickness <= 0.0:
        return vertices, faces

    min_x, min_y, min_z = vertices.min(axis=0)
    max_x, max_y, _max_z = vertices.max(axis=0)
    top_z = min_z
    bottom_z = min_z - thickness
    pad_x = max((max_x - min_x) * 0.02, 1.0)
    pad_y = max((max_y - min_y) * 0.02, 1.0)
    min_x -= pad_x
    max_x += pad_x
    min_y -= pad_y
    max_y += pad_y

    start = len(vertices)
    plate_vertices = np.asarray(
        [
            [min_x, min_y, top_z],
            [max_x, min_y, top_z],
            [max_x, max_y, top_z],
            [min_x, max_y, top_z],
            [min_x, min_y, bottom_z],
            [max_x, min_y, bottom_z],
            [max_x, max_y, bottom_z],
            [min_x, max_y, bottom_z],
        ],
        dtype=np.float64,
    )
    plate_faces = np.asarray(
        [
            [start + 0, start + 1, start + 2],
            [start + 0, start + 2, start + 3],
            [start + 4, start + 6, start + 5],
            [start + 4, start + 7, start + 6],
            [start + 0, start + 4, start + 5],
            [start + 0, start + 5, start + 1],
            [start + 1, start + 5, start + 6],
            [start + 1, start + 6, start + 2],
            [start + 2, start + 6, start + 7],
            [start + 2, start + 7, start + 3],
            [start + 3, start + 7, start + 4],
            [start + 3, start + 4, start + 0],
        ],
        dtype=np.int64,
    )
    return np.vstack([vertices, plate_vertices]), np.vstack([faces, plate_faces])


def maybe_repair(mesh) -> None:
    import trimesh

    mesh.update_faces(mesh.unique_faces())
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(mesh)
    try:
        trimesh.repair.fix_inversion(mesh)
    except Exception as exc:
        print(f"warning: repair inversion skipped: {exc}", file=sys.stderr)


def build_stl_in_memory(args: argparse.Namespace) -> Stats:
    import trimesh

    root = Path(args.root)
    out_path = Path(args.out)
    bbox = BBox(west=args.west, east=args.east, south=args.south, north=args.north)
    transforms = make_transformers()
    stats = Stats()
    all_vertices_m: list[tuple[float, float, float]] = []
    all_faces: list[tuple[int, int, int]] = []

    files = filter_files_by_size(
        find_building_files(root, args.max_files, args.file_offset),
        args.skip_files_larger_than_mb,
    )
    stats.scanned_files = len(files)

    for path in tqdm(files, desc="citygml", unit="file"):
        used_this_file = False
        try:
            for ring in iter_rings(path, transforms, args.lod):
                stats.scanned_polygons += 1
                if not bbox.contains(ring.centroid_lon, ring.centroid_lat):
                    continue
                append_ring_mesh(ring, all_vertices_m, all_faces)
                stats.used_polygons += 1
                used_this_file = True
        except Exception as exc:
            print(f"warning: skipped CityGML file: {path}: {exc}", file=sys.stderr)
            continue
        if used_this_file:
            stats.used_files += 1

    if not all_vertices_m or not all_faces:
        raise RuntimeError(
            "empty STL: no polygons matched the bbox. "
            "Check west/east/south/north, try a larger bbox, or verify the source CRS detection."
        )

    bounds = Bounds()
    for point in all_vertices_m:
        bounds.update(point)
    coordinates = coordinates_for_bounds(args, bounds)
    vertices, faces = scaled_mesh_data(all_vertices_m, all_faces, args.scale, coordinates)
    vertices, faces = add_base_plate(vertices, faces, args.base_thickness)

    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    if args.repair:
        maybe_repair(mesh)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(out_path)
    if getattr(args, "origin_out", None):
        coordinates.write(args.origin_out)

    stats.vertices = int(len(mesh.vertices))
    stats.faces = int(len(mesh.faces))
    return stats


def build_stl(args: argparse.Namespace) -> Stats:
    root = Path(args.root)
    out_path = Path(args.out)
    bbox = BBox(west=args.west, east=args.east, south=args.south, north=args.north)
    transforms = make_transformers()
    files = filter_files_by_size(
        find_building_files(root, args.max_files, args.file_offset),
        args.skip_files_larger_than_mb,
    )
    stats, bounds, face_count = scan_stl_inputs(files, bbox, transforms, args.lod)

    if stats.used_polygons == 0 or face_count == 0:
        raise RuntimeError(
            "empty STL: no polygons matched the bbox. "
            "Check west/east/south/north, try a larger bbox, or verify the source CRS detection."
        )

    coordinates = coordinates_for_bounds(args, bounds)
    base_vertices, base_faces = write_streaming_stl(
        files=files,
        out_path=out_path,
        bbox=bbox,
        transforms=transforms,
        bounds=bounds,
        face_count=face_count,
        scale=args.scale,
        base_thickness=args.base_thickness,
        lod=args.lod,
        coordinates=coordinates,
    )
    if getattr(args, "origin_out", None):
        coordinates.write(args.origin_out)
    stats.vertices += base_vertices
    stats.faces += base_faces
    return stats


def command_stl(args: argparse.Namespace) -> int:
    if args.repair or args.in_memory:
        if args.repair:
            print("warning: --repair uses in-memory mesh processing", file=sys.stderr)
        stats = build_stl_in_memory(args)
    else:
        stats = build_stl(args)
    print(f"scanned files: {stats.scanned_files}")
    print(f"used files: {stats.used_files}")
    print(f"scanned polygons: {stats.scanned_polygons}")
    print(f"used polygons: {stats.used_polygons}")
    print(f"vertices: {stats.vertices}")
    print(f"faces: {stats.faces}")
    print(f"output: {args.out}")
    return 0


def add_common_bbox_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--west", type=float, required=True)
    parser.add_argument("--east", type=float, required=True)
    parser.add_argument("--south", type=float, required=True)
    parser.add_argument("--north", type=float, required=True)


def build_parser() -> argparse.ArgumentParser:
    from subway_pipeline import command_subway_init
    from subway_pipeline import command_subway_normalize
    from subway_pipeline import command_subway_report
    from subway_pipeline import command_subway_rough_stl

    parser = argparse.ArgumentParser(description="Convert PLATEAU CityGML buildings to STL.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_download = sub.add_parser("download", help="download a CityGML ZIP file")
    p_download.add_argument("--url", required=True)
    p_download.add_argument("--out", required=True)
    p_download.set_defaults(func=command_download)

    p_unzip = sub.add_parser("unzip", help="extract a CityGML ZIP file")
    p_unzip.add_argument("--zip", required=True)
    p_unzip.add_argument("--out", required=True)
    p_unzip.set_defaults(func=command_unzip)

    p_stl = sub.add_parser("stl", help="convert building CityGML in bbox to STL")
    p_stl.add_argument("--root", required=True)
    p_stl.add_argument("--out", required=True)
    add_common_bbox_args(p_stl)
    p_stl.add_argument("--scale", "--mm-per-m", dest="scale", type=float, default=0.001,
                       help="millimetres per source metre; 1:120000 uses 0.008333333333333333")
    p_stl.add_argument("--origin-json", help="fixed ENU projected-metre/model-mm transform")
    p_stl.add_argument("--origin-out", "--export-origin", dest="origin_out",
                       help="write the exact transform used for later batches")
    p_stl.add_argument("--vertical-exaggeration", type=float, default=1.0)
    p_stl.add_argument("--z-offset", type=float, default=0.0, help="model millimetres")
    p_stl.add_argument("--base-thickness", type=float, default=0.0)
    p_stl.add_argument("--max-files", type=int, default=None)
    p_stl.add_argument("--file-offset", type=int, default=0)
    p_stl.add_argument("--lod", choices=("any", "0", "1", "2", "3"), default="any")
    p_stl.add_argument("--skip-files-larger-than-mb", type=float, default=None)
    p_stl.add_argument("--in-memory", action="store_true")
    p_stl.add_argument("--repair", action="store_true")
    p_stl.set_defaults(func=command_stl)

    p_subway_init = sub.add_parser("subway-init", help="create subway data templates")
    p_subway_init.add_argument("--out-dir", required=True)
    p_subway_init.add_argument("--sample", action=argparse.BooleanOptionalAction, default=True)
    p_subway_init.set_defaults(func=command_subway_init)

    p_subway_norm = sub.add_parser("subway-normalize", help="normalize subway source data")
    p_subway_norm.add_argument("--railway-root", default=None)
    p_subway_norm.add_argument("--gtfs-dir", action="append", default=[])
    p_subway_norm.add_argument("--lines-yaml", required=True)
    p_subway_norm.add_argument("--depth-csv", required=True)
    p_subway_norm.add_argument("--aliases-csv", default=None)
    p_subway_norm.add_argument("--out-dir", required=True)
    p_subway_norm.add_argument("--west", type=float, default=139.655)
    p_subway_norm.add_argument("--east", type=float, default=139.735)
    p_subway_norm.add_argument("--south", type=float, default=35.640)
    p_subway_norm.add_argument("--north", type=float, default=35.715)
    p_subway_norm.add_argument("--clip-mode", choices=("bbox", "all"), default="bbox")
    p_subway_norm.set_defaults(func=command_subway_normalize)

    p_subway_stl = sub.add_parser("subway-rough-stl", help="create rough subway STL")
    p_subway_stl.add_argument("--normalized-dir", required=True)
    p_subway_stl.add_argument("--origin-json", default=None)
    p_subway_stl.add_argument("--out", required=True)
    p_subway_stl.add_argument("--tube-out", default=None)
    p_subway_stl.add_argument("--shaft-out", default=None)
    p_subway_stl.add_argument("--by-line-dir", default="output/subway_rough_by_line")
    p_subway_stl.add_argument("--scale", type=float, default=0.0005)
    p_subway_stl.add_argument("--tube-radius", type=float, default=1.0)
    p_subway_stl.add_argument("--shaft-radius", type=float, default=1.5)
    p_subway_stl.add_argument("--tube-sides", type=int, default=8)
    p_subway_stl.add_argument("--shaft-sides", type=int, default=8)
    p_subway_stl.add_argument("--vertical-exaggeration", type=float, default=2.0)
    p_subway_stl.add_argument("--z-offset", type=float, default=0.0)
    p_subway_stl.add_argument("--station-shafts", action="store_true")
    p_subway_stl.add_argument("--split-by-line", action="store_true")
    p_subway_stl.add_argument("--debug-points", default=None)
    p_subway_stl.set_defaults(func=command_subway_rough_stl)

    p_subway_report = sub.add_parser("subway-report", help="create subway data quality report")
    p_subway_report.add_argument("--normalized-dir", required=True)
    p_subway_report.add_argument("--out", required=True)
    p_subway_report.set_defaults(func=command_subway_report)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "scale", 1.0) <= 0.0:
        parser.error("--scale must be greater than zero")
    if getattr(args, "vertical_exaggeration", 1.0) <= 0.0:
        parser.error("--vertical-exaggeration must be greater than zero")
    if getattr(args, "base_thickness", 0.0) < 0.0:
        parser.error("--base-thickness must be zero or greater")
    return int(args.func(args))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)

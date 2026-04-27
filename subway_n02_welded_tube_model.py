#!/usr/bin/env python3
"""Build rough subway STL from N02 railway GeoJSON with welded tube rings."""

from __future__ import annotations

import argparse
import csv
import json
import math
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests
from pyproj import Transformer


WGS84 = "EPSG:4326"
TOKYO_METERS = "EPSG:6677"
N02_URL = "https://nlftp.mlit.go.jp/ksj/gml/data/N02/N02-24/N02-24_GML.zip"


LINE_MAP = [
    ("Tokyo Metro", "Ginza Line", "東京地下鉄", "3号線銀座線", 12.0),
    ("Tokyo Metro", "Marunouchi Line", "東京地下鉄", "4号線丸ノ内線", 18.0),
    ("Tokyo Metro", "Marunouchi Line", "東京地下鉄", "4号線丸ノ内線分岐線", 18.0),
    ("Tokyo Metro", "Hibiya Line", "東京地下鉄", "2号線日比谷線", 18.0),
    ("Tokyo Metro", "Tozai Line", "東京地下鉄", "5号線東西線", 20.0),
    ("Tokyo Metro", "Chiyoda Line", "東京地下鉄", "9号線千代田線", 23.0),
    ("Tokyo Metro", "Yurakucho Line", "東京地下鉄", "8号線有楽町線", 25.0),
    ("Tokyo Metro", "Hanzomon Line", "東京地下鉄", "11号線半蔵門線", 28.0),
    ("Tokyo Metro", "Namboku Line", "東京地下鉄", "7号線南北線", 28.0),
    ("Tokyo Metro", "Fukutoshin Line", "東京地下鉄", "13号線副都心線", 30.0),
    ("Toei Subway", "Asakusa Line", "東京都", "1号線浅草線", 20.0),
    ("Toei Subway", "Mita Line", "東京都", "6号線三田線", 24.0),
    ("Toei Subway", "Shinjuku Line", "東京都", "10号線新宿線", 25.0),
    ("Toei Subway", "Oedo Line", "東京都", "12号線大江戸線", 35.0),
]


@dataclass
class Origin:
    origin_x: float
    origin_y: float
    origin_z: float
    scale: float


def ensure_n02_zip(path: Path) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(N02_URL, timeout=120)
    response.raise_for_status()
    path.write_bytes(response.content)
    return path


def load_geojson_from_zip(path: Path, member: str) -> dict[str, Any]:
    with zipfile.ZipFile(path, "r") as zf:
        return json.loads(zf.read(member).decode("utf-8"))


def line_lookup() -> dict[tuple[str, str], tuple[str, str, float]]:
    return {(operator_jp, line_jp): (operator, line, depth) for operator, line, operator_jp, line_jp, depth in LINE_MAP}


def bbox_contains_coord(coord: list[float], args: argparse.Namespace) -> bool:
    lon = float(coord[0])
    lat = float(coord[1])
    return args.west <= lon <= args.east and args.south <= lat <= args.north


def line_intersects_bbox(coords: list[list[float]], args: argparse.Namespace) -> bool:
    if args.clip_mode == "all":
        return True
    return any(bbox_contains_coord(coord, args) for coord in coords)


def point_for_station_geometry(geometry: dict[str, Any]) -> list[float] | None:
    coords = geometry.get("coordinates") or []
    if geometry.get("type") == "Point":
        return coords
    if geometry.get("type") == "LineString" and coords:
        mid = coords[len(coords) // 2]
        return mid
    return None


def load_n02(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    path = ensure_n02_zip(Path(args.n02_zip))
    lookup = line_lookup()
    line_data = load_geojson_from_zip(path, "UTF-8/N02-24_RailroadSection.geojson")
    station_data = load_geojson_from_zip(path, "UTF-8/N02-24_Station.geojson")
    lines = []
    stations = []

    for feature in line_data.get("features", []):
        props = feature.get("properties", {})
        key = (props.get("N02_004", ""), props.get("N02_003", ""))
        canon = lookup.get(key)
        if canon is None:
            continue
        geometry = feature.get("geometry") or {}
        if geometry.get("type") != "LineString":
            continue
        coords = geometry.get("coordinates") or []
        if len(coords) < 2 or not line_intersects_bbox(coords, args):
            continue
        lines.append(
            {
                "canonical_operator": canon[0],
                "canonical_line": canon[1],
                "default_depth_m": canon[2],
                "coords": coords,
                "source_line_name": key[1],
                "source_operator": key[0],
            }
        )

    for feature in station_data.get("features", []):
        props = feature.get("properties", {})
        key = (props.get("N02_004", ""), props.get("N02_003", ""))
        canon = lookup.get(key)
        if canon is None:
            continue
        coord = point_for_station_geometry(feature.get("geometry") or {})
        if coord is None:
            continue
        if args.clip_mode != "all" and not bbox_contains_coord(coord, args):
            continue
        stations.append(
            {
                "canonical_operator": canon[0],
                "canonical_line": canon[1],
                "canonical_station": props.get("N02_005", ""),
                "station_id": props.get("N02_005c", "") or props.get("N02_005g", ""),
                "lon": float(coord[0]),
                "lat": float(coord[1]),
                "surface_z_m": 0.0,
                "depth_m": canon[2],
                "station_z_m": -canon[2],
                "depth_type": "estimated",
                "confidence": "estimated",
                "source_label": "N02-24",
                "source_url": N02_URL,
                "note": "N02 station position with route default depth.",
            }
        )
    return lines, dedupe_stations(stations)


def dedupe_stations(stations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen = set()
    out = []
    for row in sorted(stations, key=lambda r: (r["canonical_operator"], r["canonical_line"], r["canonical_station"], r["station_id"])):
        key = (row["canonical_operator"], row["canonical_line"], row["canonical_station"], round(row["lon"], 6), round(row["lat"], 6))
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def add_epsg6677(lines: list[dict[str, Any]], stations: list[dict[str, Any]]) -> None:
    transformer = Transformer.from_crs(WGS84, TOKYO_METERS, always_xy=True)
    for line in lines:
        xy = []
        for lon, lat in line["coords"]:
            x, y = transformer.transform(lon, lat)
            xy.append([x, y])
        line["xy"] = xy
    for row in stations:
        x, y = transformer.transform(row["lon"], row["lat"])
        row["x_epsg6677"] = x
        row["y_epsg6677"] = y


def make_origin(lines: list[dict[str, Any]], stations: list[dict[str, Any]], args: argparse.Namespace) -> Origin:
    if args.origin_json:
        data = json.loads(Path(args.origin_json).read_text(encoding="utf-8"))
        return Origin(float(data["origin_x"]), float(data["origin_y"]), float(data.get("origin_z", 0.0)), float(data.get("scale", args.scale)))
    xs = []
    ys = []
    for line in lines:
        xs.extend(point[0] for point in line["xy"])
        ys.extend(point[1] for point in line["xy"])
    for row in stations:
        xs.append(float(row["x_epsg6677"]))
        ys.append(float(row["y_epsg6677"]))
    return Origin((min(xs) + max(xs)) * 0.5, (min(ys) + max(ys)) * 0.5, 0.0, args.scale)


def model_xy(x: float, y: float, z_m: float, origin: Origin, args: argparse.Namespace) -> list[float]:
    return [
        (x - origin.origin_x) * origin.scale,
        (y - origin.origin_y) * origin.scale,
        (z_m - origin.origin_z) * origin.scale * args.vertical_exaggeration + args.z_offset,
    ]


def line_model_points(line: dict[str, Any], origin: Origin, args: argparse.Namespace) -> list[list[float]]:
    z_m = -float(line["default_depth_m"])
    return [model_xy(x, y, z_m, origin, args) for x, y in simplify_xy(line["xy"], args.min_point_spacing_m)]


def xy_distance(a: list[float], b: list[float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def segment_projection_xy(point: list[float], a: list[float], b: list[float]) -> tuple[float, float, float]:
    abx = b[0] - a[0]
    aby = b[1] - a[1]
    apx = point[0] - a[0]
    apy = point[1] - a[1]
    denom = abx * abx + aby * aby
    if denom <= 0.0:
        return xy_distance(point, a), a[0], a[1]
    t = max(0.0, min(1.0, (apx * abx + apy * aby) / denom))
    x = a[0] + abx * t
    y = a[1] + aby * t
    return math.hypot(point[0] - x, point[1] - y), x, y


def nearest_line_yaw(point_xy: list[float], lines: list[dict[str, Any]], fallback_yaw: float = 0.0) -> float:
    best_dist = math.inf
    best_yaw = fallback_yaw
    for line in lines:
        xy = line.get("xy", [])
        for idx in range(len(xy) - 1):
            a = xy[idx]
            b = xy[idx + 1]
            dist, _x, _y = segment_projection_xy(point_xy, a, b)
            if dist >= best_dist:
                continue
            dx = b[0] - a[0]
            dy = b[1] - a[1]
            if dx == 0.0 and dy == 0.0:
                continue
            best_dist = dist
            best_yaw = math.atan2(dy, dx)
    return best_yaw


def stitch_line_parts(lines: list[dict[str, Any]], tolerance_m: float) -> list[dict[str, Any]]:
    if tolerance_m <= 0.0:
        return lines
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for line in lines:
        grouped[(line["canonical_operator"], line["canonical_line"])].append(line)

    stitched = []
    for (_operator, _line), parts in grouped.items():
        remaining = [dict(part) for part in parts if len(part.get("xy", [])) >= 2]
        while remaining:
            current = remaining.pop(0)
            path = list(current["xy"])
            changed = True
            while changed:
                changed = False
                for idx, candidate in enumerate(remaining):
                    cand = list(candidate["xy"])
                    joins = [
                        (xy_distance(path[-1], cand[0]), "append"),
                        (xy_distance(path[-1], cand[-1]), "append_reverse"),
                        (xy_distance(path[0], cand[-1]), "prepend"),
                        (xy_distance(path[0], cand[0]), "prepend_reverse"),
                    ]
                    dist, mode = min(joins, key=lambda item: item[0])
                    if dist > tolerance_m:
                        continue
                    if mode == "append":
                        path.extend(cand[1:])
                    elif mode == "append_reverse":
                        path.extend(list(reversed(cand[:-1])))
                    elif mode == "prepend":
                        path = cand[:-1] + path
                    else:
                        path = list(reversed(cand[1:])) + path
                    remaining.pop(idx)
                    changed = True
                    break
            merged = dict(current)
            merged["xy"] = path
            stitched.append(merged)
    return stitched


def simplify_xy(points: list[list[float]], min_spacing: float) -> list[list[float]]:
    if len(points) <= 2 or min_spacing <= 0.0:
        return points
    out = [points[0]]
    last = points[0]
    for point in points[1:-1]:
        if math.hypot(point[0] - last[0], point[1] - last[1]) >= min_spacing:
            out.append(point)
            last = point
    out.append(points[-1])
    return out


def make_cylinder(p1: list[float], p2: list[float], radius: float, sides: int) -> Any | None:
    import numpy as np
    import trimesh

    if np.linalg.norm(np.asarray(p2) - np.asarray(p1)) <= 0.0:
        return None
    return trimesh.creation.cylinder(radius=radius, sections=sides, segment=[p1, p2])


def make_box_at(center: list[float], size: float) -> Any:
    import trimesh

    mesh = trimesh.creation.box(extents=[size, size, size])
    mesh.apply_translation(center)
    return mesh


def make_oriented_box_at(center: list[float], yaw: float, length: float, width: float, height: float) -> Any:
    import numpy as np
    import trimesh

    mesh = trimesh.creation.box(extents=[length, width, height])
    c = math.cos(yaw)
    s = math.sin(yaw)
    transform = np.array(
        [
            [c, -s, 0.0, center[0]],
            [s, c, 0.0, center[1]],
            [0.0, 0.0, 1.0, center[2]],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    mesh.apply_transform(transform)
    return mesh


def vector_sub(a: list[float], b: list[float]) -> list[float]:
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def vector_add(a: list[float], b: list[float]) -> list[float]:
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def vector_scale(a: list[float], scale: float) -> list[float]:
    return [a[0] * scale, a[1] * scale, a[2] * scale]


def dot(a: list[float], b: list[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: list[float], b: list[float]) -> list[float]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def norm(a: list[float]) -> float:
    return math.sqrt(dot(a, a))


def normalize(a: list[float]) -> list[float]:
    length = norm(a)
    if length <= 0.0:
        return [0.0, 0.0, 0.0]
    return [a[0] / length, a[1] / length, a[2] / length]


def tangent_at(points: list[list[float]], index: int) -> list[float]:
    if index == 0:
        return normalize(vector_sub(points[1], points[0]))
    if index == len(points) - 1:
        return normalize(vector_sub(points[-1], points[-2]))
    before = normalize(vector_sub(points[index], points[index - 1]))
    after = normalize(vector_sub(points[index + 1], points[index]))
    mixed = vector_add(before, after)
    if norm(mixed) <= 0.0:
        return after
    return normalize(mixed)


def frame_for_tangent(tangent: list[float]) -> tuple[list[float], list[float]]:
    up = [0.0, 0.0, 1.0]
    if abs(dot(tangent, up)) > 0.95:
        up = [1.0, 0.0, 0.0]
    normal = normalize(cross(tangent, up))
    binormal = normalize(cross(tangent, normal))
    return normal, binormal


def make_welded_tube(points: list[list[float]], radius: float, sides: int) -> Any | None:
    if len(points) < 2:
        return None
    import numpy as np
    import trimesh

    clean_points = [points[0]]
    for point in points[1:]:
        if norm(vector_sub(point, clean_points[-1])) > 0.0:
            clean_points.append(point)
    if len(clean_points) < 2:
        return None

    vertices = []
    faces = []
    for index, point in enumerate(clean_points):
        tangent = tangent_at(clean_points, index)
        normal, binormal = frame_for_tangent(tangent)
        for side in range(sides):
            angle = 2.0 * math.pi * side / sides
            offset = vector_add(
                vector_scale(normal, math.cos(angle) * radius),
                vector_scale(binormal, math.sin(angle) * radius),
            )
            vertices.append(vector_add(point, offset))

    for ring in range(len(clean_points) - 1):
        base_a = ring * sides
        base_b = (ring + 1) * sides
        for side in range(sides):
            a0 = base_a + side
            a1 = base_a + ((side + 1) % sides)
            b0 = base_b + side
            b1 = base_b + ((side + 1) % sides)
            faces.append([a0, b0, b1])
            faces.append([a0, b1, a1])

    start_center = len(vertices)
    vertices.append(clean_points[0])
    end_center = len(vertices)
    vertices.append(clean_points[-1])
    for side in range(sides):
        next_side = (side + 1) % sides
        faces.append([start_center, next_side, side])
        end_base = (len(clean_points) - 1) * sides
        faces.append([end_center, end_base + side, end_base + next_side])

    return trimesh.Trimesh(vertices=np.asarray(vertices), faces=np.asarray(faces), process=False)


def concatenate(meshes: list[Any]) -> Any:
    if not meshes:
        return None
    import trimesh

    return trimesh.util.concatenate(meshes)


def export_mesh(mesh: Any, path: Path) -> None:
    if mesh is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(path)


def write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def write_geojson(path: Path, features: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}, ensure_ascii=False, indent=2), encoding="utf-8")


def write_outputs(lines: list[dict[str, Any]], stations: list[dict[str, Any]], args: argparse.Namespace) -> None:
    out_dir = Path(args.normalized_out_dir)
    station_cols = [
        "canonical_operator",
        "canonical_line",
        "canonical_station",
        "station_id",
        "lon",
        "lat",
        "x_epsg6677",
        "y_epsg6677",
        "surface_z_m",
        "depth_m",
        "station_z_m",
        "depth_type",
        "confidence",
        "source_label",
        "source_url",
        "note",
    ]
    write_csv(out_dir / "n02_subway_stations_epsg6677.csv", stations, station_cols)
    write_geojson(
        out_dir / "n02_subway_stations.geojson",
        [
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [row["lon"], row["lat"]]},
                "properties": {col: row.get(col, "") for col in station_cols},
            }
            for row in stations
        ],
    )
    write_geojson(
        out_dir / "n02_subway_lines.geojson",
        [
            {
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": line["coords"]},
                "properties": {
                    "canonical_operator": line["canonical_operator"],
                    "canonical_line": line["canonical_line"],
                    "source_line_name": line["source_line_name"],
                    "source_operator": line["source_operator"],
                    "default_depth_m": line["default_depth_m"],
                },
            }
            for line in lines
        ],
    )


def normalize_name(value: str) -> str:
    out = "".join(ch.lower() if ch.isalnum() else "_" for ch in value)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def build(args: argparse.Namespace) -> int:
    import trimesh  # noqa: F401

    lines, stations = load_n02(args)
    add_epsg6677(lines, stations)
    origin = make_origin(lines, stations, args)
    write_outputs(lines, stations, args)

    tube_radius = args.tube_radius_m * origin.scale
    shaft_radius = args.shaft_radius_m * origin.scale
    tube_meshes = []
    shaft_meshes = []
    junction_meshes = []
    by_line: dict[tuple[str, str], list[Any]] = defaultdict(list)

    mesh_lines = stitch_line_parts(lines, args.stitch_tolerance_m)

    for line in mesh_lines:
        points = line_model_points(line, origin, args)
        mesh = make_welded_tube(points, tube_radius, args.tube_sides)
        if mesh is None:
            continue
        tube_meshes.append(mesh)
        by_line[(line["canonical_operator"], line["canonical_line"])].append(mesh)

    if args.station_shafts:
        for row in stations:
            p1 = model_xy(float(row["x_epsg6677"]), float(row["y_epsg6677"]), float(row["station_z_m"]), origin, args)
            p2 = model_xy(float(row["x_epsg6677"]), float(row["y_epsg6677"]), float(row["surface_z_m"]), origin, args)
            mesh = make_cylinder(p1, p2, shaft_radius, args.shaft_sides)
            if mesh is not None:
                shaft_meshes.append(mesh)
            if args.junction_boxes:
                yaw = nearest_line_yaw([float(row["x_epsg6677"]), float(row["y_epsg6677"])], mesh_lines)
                if args.junction_box_size_m > 0.0:
                    box_size = args.junction_box_size_m * origin.scale
                    junction_meshes.append(make_oriented_box_at(p1, yaw, box_size, box_size, box_size))
                else:
                    junction_meshes.append(
                        make_oriented_box_at(
                            p1,
                            yaw,
                            args.junction_box_length_m * origin.scale,
                            args.junction_box_width_m * origin.scale,
                            args.junction_box_height_m * origin.scale,
                        )
                    )

    tube_mesh = concatenate(tube_meshes)
    shaft_mesh = concatenate(shaft_meshes)
    junction_mesh = concatenate(junction_meshes)
    combined = concatenate([mesh for mesh in [tube_mesh, shaft_mesh, junction_mesh] if mesh is not None])
    export_mesh(tube_mesh, Path(args.tube_out))
    export_mesh(shaft_mesh, Path(args.shaft_out))
    if args.junction_box_out:
        export_mesh(junction_mesh, Path(args.junction_box_out))
    export_mesh(combined, Path(args.out))
    if args.split_by_line:
        out_dir = Path(args.by_line_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for (operator, line), meshes in by_line.items():
            export_mesh(concatenate(meshes), out_dir / f"{normalize_name(operator + '_' + line)}.stl")

    counts: dict[str, int] = defaultdict(int)
    for station in stations:
        counts[f"{station['canonical_operator']} {station['canonical_line']}"] += 1
    print(f"n02 lines: {len(lines)}")
    print(f"stitched tube paths: {len(mesh_lines)}")
    print(f"n02 stations: {len(stations)}")
    print(f"tube meshes: {len(tube_meshes)}")
    print(f"shafts: {len(shaft_meshes)}")
    print(f"junction boxes: {len(junction_meshes)}")
    for name in sorted(counts):
        print(f"{name}: {counts[name]}")
    print(f"tube_radius_model: {tube_radius}")
    print(f"output: {args.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build welded subway rough STL from N02 line geometry.")
    parser.add_argument("--n02-zip", default="data/n02/N02-24_GML.zip")
    parser.add_argument("--normalized-out-dir", default="normalized_subway_n02")
    parser.add_argument("--origin-json", default=None)
    parser.add_argument("--out", default="output/subway_n02_rough_with_shafts.stl")
    parser.add_argument("--tube-out", default="output/subway_n02_rough_all_lines.stl")
    parser.add_argument("--shaft-out", default="output/subway_n02_station_shafts.stl")
    parser.add_argument("--junction-box-out", default=None)
    parser.add_argument("--by-line-dir", default="output/subway_n02_by_line")
    parser.add_argument("--west", type=float, default=139.655)
    parser.add_argument("--east", type=float, default=139.735)
    parser.add_argument("--south", type=float, default=35.640)
    parser.add_argument("--north", type=float, default=35.715)
    parser.add_argument("--clip-mode", choices=("bbox", "all"), default="bbox")
    parser.add_argument("--scale", type=float, default=0.0005)
    parser.add_argument("--tube-radius-m", type=float, default=10.0)
    parser.add_argument("--shaft-radius-m", type=float, default=16.0)
    parser.add_argument("--junction-box-size-m", type=float, default=80.0)
    parser.add_argument("--junction-box-length-m", type=float, default=180.0)
    parser.add_argument("--junction-box-width-m", type=float, default=140.0)
    parser.add_argument("--junction-box-height-m", type=float, default=40.0)
    parser.add_argument("--tube-sides", type=int, default=8)
    parser.add_argument("--shaft-sides", type=int, default=8)
    parser.add_argument("--min-point-spacing-m", type=float, default=30.0)
    parser.add_argument("--stitch-tolerance-m", type=float, default=5.0)
    parser.add_argument("--vertical-exaggeration", type=float, default=2.0)
    parser.add_argument("--z-offset", type=float, default=0.0)
    parser.add_argument("--station-shafts", action="store_true")
    parser.add_argument("--junction-boxes", action="store_true")
    parser.add_argument("--split-by-line", action="store_true")
    return parser


def main() -> int:
    return build(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build a Yamanote Line overlay STL from N02 data."""

from __future__ import annotations

import argparse
import json
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests
from pyproj import Transformer


WGS84 = "EPSG:4326"
TOKYO_METERS = "EPSG:6677"
N02_URL = "https://nlftp.mlit.go.jp/ksj/gml/data/N02/N02-24/N02-24_GML.zip"


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


def load_origin(path: Path, scale: float) -> Origin:
    data = json.loads(path.read_text(encoding="utf-8"))
    return Origin(
        origin_x=float(data.get("origin_x", 0.0)),
        origin_y=float(data.get("origin_y", 0.0)),
        origin_z=float(data.get("origin_z", 0.0)),
        scale=float(data.get("scale", scale)),
    )


def in_bbox(coord: list[float], args: argparse.Namespace) -> bool:
    lon = float(coord[0])
    lat = float(coord[1])
    return args.west <= lon <= args.east and args.south <= lat <= args.north


def load_yamanote(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    path = ensure_n02_zip(Path(args.n02_zip))
    line_data = load_geojson_from_zip(path, "UTF-8/N02-24_RailroadSection.geojson")
    station_data = load_geojson_from_zip(path, "UTF-8/N02-24_Station.geojson")
    lines = []
    stations = []
    for feature in line_data.get("features", []):
        props = feature.get("properties", {})
        if props.get("N02_003") != "山手線" or props.get("N02_004") != "東日本旅客鉄道":
            continue
        coords = (feature.get("geometry") or {}).get("coordinates") or []
        if len(coords) < 2:
            continue
        if args.clip_mode != "all" and not any(in_bbox(coord, args) for coord in coords):
            continue
        lines.append({"coords": coords, "line": "Yamanote Line", "operator": "JR East"})
    for feature in station_data.get("features", []):
        props = feature.get("properties", {})
        if props.get("N02_003") != "山手線" or props.get("N02_004") != "東日本旅客鉄道":
            continue
        coords = (feature.get("geometry") or {}).get("coordinates") or []
        if not coords:
            continue
        if (feature.get("geometry") or {}).get("type") == "LineString":
            coord = coords[len(coords) // 2]
        else:
            coord = coords
        if args.clip_mode != "all" and not in_bbox(coord, args):
            continue
        stations.append(
            {
                "station": props.get("N02_005", ""),
                "station_id": props.get("N02_005c", "") or props.get("N02_005g", ""),
                "lon": float(coord[0]),
                "lat": float(coord[1]),
            }
        )
    return lines, dedupe_stations(stations)


def dedupe_stations(stations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen = set()
    out = []
    for row in sorted(stations, key=lambda item: (item["station"], item["station_id"])):
        key = (row["station"], round(row["lon"], 6), round(row["lat"], 6))
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def add_xy(lines: list[dict[str, Any]], stations: list[dict[str, Any]]) -> None:
    transformer = Transformer.from_crs(WGS84, TOKYO_METERS, always_xy=True)
    for line in lines:
        line["xy"] = [list(transformer.transform(lon, lat)) for lon, lat in line["coords"]]
    for row in stations:
        row["xy"] = list(transformer.transform(row["lon"], row["lat"]))


def model_point(x: float, y: float, z_m: float, origin: Origin, args: argparse.Namespace) -> list[float]:
    return [
        (x - origin.origin_x) * origin.scale,
        (y - origin.origin_y) * origin.scale,
        (z_m - origin.origin_z) * origin.scale * args.vertical_exaggeration + args.z_offset,
    ]


def dist_xy(a: list[float], b: list[float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def simplify(points: list[list[float]], min_spacing: float) -> list[list[float]]:
    if len(points) <= 2:
        return points
    out = [points[0]]
    last = points[0]
    for point in points[1:-1]:
        if dist_xy(point, last) >= min_spacing:
            out.append(point)
            last = point
    out.append(points[-1])
    return out


def stitch_lines(lines: list[dict[str, Any]], tolerance_m: float) -> list[list[list[float]]]:
    remaining = [line["xy"] for line in lines if len(line.get("xy", [])) >= 2]
    paths = []
    while remaining:
        path = list(remaining.pop(0))
        changed = True
        while changed:
            changed = False
            for idx, cand in enumerate(remaining):
                joins = [
                    (dist_xy(path[-1], cand[0]), "append"),
                    (dist_xy(path[-1], cand[-1]), "append_reverse"),
                    (dist_xy(path[0], cand[-1]), "prepend"),
                    (dist_xy(path[0], cand[0]), "prepend_reverse"),
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
        paths.append(path)
    return paths


def dot(a: list[float], b: list[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: list[float], b: list[float]) -> list[float]:
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def norm(a: list[float]) -> float:
    return math.sqrt(dot(a, a))


def normalize(a: list[float]) -> list[float]:
    length = norm(a)
    if length <= 0.0:
        return [0.0, 0.0, 0.0]
    return [a[0] / length, a[1] / length, a[2] / length]


def sub(a: list[float], b: list[float]) -> list[float]:
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def add(a: list[float], b: list[float]) -> list[float]:
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def scale_vec(a: list[float], scale: float) -> list[float]:
    return [a[0] * scale, a[1] * scale, a[2] * scale]


def tangent(points: list[list[float]], idx: int) -> list[float]:
    if idx == 0:
        return normalize(sub(points[1], points[0]))
    if idx == len(points) - 1:
        return normalize(sub(points[-1], points[-2]))
    return normalize(add(normalize(sub(points[idx], points[idx - 1])), normalize(sub(points[idx + 1], points[idx]))))


def frame(tan: list[float]) -> tuple[list[float], list[float]]:
    up = [0.0, 0.0, 1.0]
    n = normalize(cross(tan, up))
    b = normalize(cross(tan, n))
    return n, b


def make_tube(points: list[list[float]], radius: float, sides: int) -> Any:
    import numpy as np
    import trimesh

    vertices = []
    faces = []
    for idx, point in enumerate(points):
        tan = tangent(points, idx)
        n, b = frame(tan)
        for side in range(sides):
            angle = 2.0 * math.pi * side / sides
            offset = add(scale_vec(n, math.cos(angle) * radius), scale_vec(b, math.sin(angle) * radius))
            vertices.append(add(point, offset))
    for ring in range(len(points) - 1):
        base_a = ring * sides
        base_b = (ring + 1) * sides
        for side in range(sides):
            a0 = base_a + side
            a1 = base_a + ((side + 1) % sides)
            b0 = base_b + side
            b1 = base_b + ((side + 1) % sides)
            faces.append([a0, b0, b1])
            faces.append([a0, b1, a1])
    return trimesh.Trimesh(vertices=np.asarray(vertices), faces=np.asarray(faces), process=False)


def make_station_markers(stations: list[dict[str, Any]], origin: Origin, args: argparse.Namespace) -> Any:
    import trimesh

    meshes = []
    radius = args.station_radius_m * origin.scale
    for row in stations:
        x, y = row["xy"]
        center = model_point(x, y, args.z_m, origin, args)
        mesh = trimesh.creation.icosphere(subdivisions=1, radius=radius)
        mesh.apply_translation(center)
        meshes.append(mesh)
    return trimesh.util.concatenate(meshes) if meshes else None


def concatenate(meshes: list[Any]) -> Any:
    import trimesh

    meshes = [mesh for mesh in meshes if mesh is not None]
    if not meshes:
        return None
    return trimesh.util.concatenate(meshes)


def export(mesh: Any, path: Path) -> None:
    if mesh is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(path)


def build(args: argparse.Namespace) -> int:
    lines, stations = load_yamanote(args)
    add_xy(lines, stations)
    origin = load_origin(Path(args.origin_json), args.scale)
    paths = stitch_lines(lines, args.stitch_tolerance_m)
    tube_meshes = []
    for path in paths:
        points = [model_point(x, y, args.z_m, origin, args) for x, y in simplify(path, args.min_point_spacing_m)]
        if len(points) >= 2:
            tube_meshes.append(make_tube(points, args.tube_radius_m * origin.scale, args.tube_sides))
    tube_mesh = concatenate(tube_meshes)
    station_mesh = make_station_markers(stations, origin, args) if args.station_markers else None
    yamanote_mesh = concatenate([tube_mesh, station_mesh])
    export(tube_mesh, Path(args.out))
    export(station_mesh, Path(args.station_out))
    export(yamanote_mesh, Path(args.with_stations_out))
    if args.subway_stl:
        import trimesh

        subway = trimesh.load_mesh(args.subway_stl, process=False)
        export(concatenate([subway, yamanote_mesh]), Path(args.combined_out))
    print(f"yamanote line parts: {len(lines)}")
    print(f"stitched paths: {len(paths)}")
    print(f"stations: {len(stations)}")
    print(f"output: {args.out}")
    if args.subway_stl:
        print(f"combined: {args.combined_out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build Yamanote Line overlay STL from N02.")
    parser.add_argument("--n02-zip", default="data/n02/N02-24_GML.zip")
    parser.add_argument("--origin-json", default="origin_epsg6677_zero.json")
    parser.add_argument("--subway-stl", default="output/subway_n02_tokyo23_welded_z20_r3_with_shafts_boxes_thin.stl")
    parser.add_argument("--out", default="output/yamanote_n02_tokyo23_overlay.stl")
    parser.add_argument("--station-out", default="output/yamanote_n02_tokyo23_stations.stl")
    parser.add_argument("--with-stations-out", default="output/yamanote_n02_tokyo23_overlay_with_stations.stl")
    parser.add_argument("--combined-out", default="output/subway_yamanote_n02_tokyo23_combined.stl")
    parser.add_argument("--west", type=float, default=139.56)
    parser.add_argument("--east", type=float, default=139.93)
    parser.add_argument("--south", type=float, default=35.52)
    parser.add_argument("--north", type=float, default=35.82)
    parser.add_argument("--clip-mode", choices=("bbox", "all"), default="bbox")
    parser.add_argument("--scale", type=float, default=0.0005)
    parser.add_argument("--z-m", type=float, default=18.0)
    parser.add_argument("--vertical-exaggeration", type=float, default=20.0)
    parser.add_argument("--z-offset", type=float, default=0.0)
    parser.add_argument("--tube-radius-m", type=float, default=24.0)
    parser.add_argument("--station-radius-m", type=float, default=45.0)
    parser.add_argument("--tube-sides", type=int, default=16)
    parser.add_argument("--min-point-spacing-m", type=float, default=50.0)
    parser.add_argument("--stitch-tolerance-m", type=float, default=8.0)
    parser.add_argument("--station-markers", action="store_true")
    return parser


def main() -> int:
    return build(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

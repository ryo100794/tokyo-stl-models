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
JR_EAST = "\u6771\u65e5\u672c\u65c5\u5ba2\u9244\u9053"
ROUTE_YAMANOTE = "\u5c71\u624b\u7dda"
ROUTE_TOHOKU = "\u6771\u5317\u7dda"
ROUTE_TOKAIDO = "\u6771\u6d77\u9053\u7dda"
YAMANOTE_STATIONS = [
    "\u6771\u4eac",
    "\u795e\u7530",
    "\u79cb\u8449\u539f",
    "\u5fa1\u5f92\u753a",
    "\u4e0a\u91ce",
    "\u9daf\u8c37",
    "\u65e5\u66ae\u91cc",
    "\u897f\u65e5\u66ae\u91cc",
    "\u7530\u7aef",
    "\u99d2\u8fbc",
    "\u5de3\u9d28",
    "\u5927\u585a",
    "\u6c60\u888b",
    "\u76ee\u767d",
    "\u9ad8\u7530\u99ac\u5834",
    "\u65b0\u5927\u4e45\u4fdd",
    "\u65b0\u5bbf",
    "\u4ee3\u3005\u6728",
    "\u539f\u5bbf",
    "\u6e0b\u8c37",
    "\u6075\u6bd4\u5bff",
    "\u76ee\u9ed2",
    "\u4e94\u53cd\u7530",
    "\u5927\u5d0e",
    "\u54c1\u5ddd",
    "\u9ad8\u8f2a\u30b2\u30fc\u30c8\u30a6\u30a7\u30a4",
    "\u7530\u753a",
    "\u6d5c\u677e\u753a",
    "\u65b0\u6a4b",
    "\u6709\u697d\u753a",
]
YAMANOTE_STATION_SET = set(YAMANOTE_STATIONS)
YAMANOTE_ROUTES = {ROUTE_YAMANOTE, ROUTE_TOHOKU, ROUTE_TOKAIDO}
SOUTHEAST_STATIONS = {
    "\u6771\u4eac",
    "\u54c1\u5ddd",
    "\u9ad8\u8f2a\u30b2\u30fc\u30c8\u30a6\u30a7\u30a4",
    "\u7530\u753a",
    "\u6d5c\u677e\u753a",
    "\u65b0\u6a4b",
    "\u6709\u697d\u753a",
}
NORTHEAST_STATIONS = set(YAMANOTE_STATIONS[:9])


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


def coord_for_feature(feature: dict[str, Any]) -> list[float] | None:
    coords = (feature.get("geometry") or {}).get("coordinates") or []
    if not coords:
        return None
    if (feature.get("geometry") or {}).get("type") == "LineString":
        return coords[len(coords) // 2]
    return coords


def station_route_preference(station: str) -> list[str]:
    if station in SOUTHEAST_STATIONS:
        return [ROUTE_TOKAIDO, ROUTE_TOHOKU, ROUTE_YAMANOTE]
    if station in NORTHEAST_STATIONS:
        return [ROUTE_TOHOKU, ROUTE_YAMANOTE, ROUTE_TOKAIDO]
    return [ROUTE_YAMANOTE, ROUTE_TOHOKU, ROUTE_TOKAIDO]


def average_candidate(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    lon = sum(row["lon"] for row in candidates) / len(candidates)
    lat = sum(row["lat"] for row in candidates) / len(candidates)
    first = candidates[0]
    return {
        "station": first["station"],
        "station_id": first["station_id"],
        "lon": lon,
        "lat": lat,
        "source_line": first["source_line"],
        "source_count": len(candidates),
    }


def choose_station(station: str, candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not candidates:
        return None
    for route in station_route_preference(station):
        rows = [row for row in candidates if row["source_line"] == route]
        if rows:
            return average_candidate(rows)
    return average_candidate(candidates)


def load_yamanote(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    path = ensure_n02_zip(Path(args.n02_zip))
    station_data = load_geojson_from_zip(path, "UTF-8/N02-24_Station.geojson")
    by_station: dict[str, list[dict[str, Any]]] = {name: [] for name in YAMANOTE_STATIONS}
    for feature in station_data.get("features", []):
        props = feature.get("properties", {})
        station = props.get("N02_005", "")
        route = props.get("N02_003", "")
        if props.get("N02_004") != JR_EAST:
            continue
        if station not in YAMANOTE_STATION_SET or route not in YAMANOTE_ROUTES:
            continue
        coord = coord_for_feature(feature)
        if coord is None:
            continue
        if args.clip_mode != "all" and not in_bbox(coord, args):
            continue
        by_station[station].append(
            {
                "station": station,
                "station_id": props.get("N02_005c", "") or props.get("N02_005g", ""),
                "lon": float(coord[0]),
                "lat": float(coord[1]),
                "source_line": route,
            }
        )
    stations = []
    missing = []
    for station in YAMANOTE_STATIONS:
        row = choose_station(station, by_station.get(station, []))
        if row is None:
            missing.append(station)
            continue
        stations.append(row)
    if missing:
        missing_text = ", ".join(missing)
        raise RuntimeError(f"missing Yamanote stations in N02 data: {missing_text}")
    line_coords = [[row["lon"], row["lat"]] for row in stations]
    line_coords.append([stations[0]["lon"], stations[0]["lat"]])
    lines = [{"coords": line_coords, "line": "Yamanote Line", "operator": "JR East"}]
    return lines, stations


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

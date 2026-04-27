#!/usr/bin/env python3
"""Build rough geology layer STL surfaces for Tokyo 23 wards."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import trimesh
from pyproj import Transformer


WGS84 = "EPSG:4326"
TOKYO_METERS = "EPSG:6677"


@dataclass
class Origin:
    origin_x: float
    origin_y: float
    origin_z: float
    scale: float


DEFAULT_BOUNDARIES = [
    ("surface", 0.0),
    ("base_fill", -2.0),
    ("base_alluvium", -18.0),
    ("base_tokyo_upper", -42.0),
    ("base_tokyo_lower", -78.0),
    ("base_kazusa_top", -125.0),
]


LAYER_PAIRS = [
    ("fill", "surface", "base_fill"),
    ("alluvium", "base_fill", "base_alluvium"),
    ("tokyo_upper", "base_alluvium", "base_tokyo_upper"),
    ("tokyo_lower", "base_tokyo_upper", "base_tokyo_lower"),
    ("kazusa_top", "base_tokyo_lower", "base_kazusa_top"),
]


def load_origin(path: Path, scale: float) -> Origin:
    data = json.loads(path.read_text(encoding="utf-8"))
    return Origin(
        origin_x=float(data.get("origin_x", 0.0)),
        origin_y=float(data.get("origin_y", 0.0)),
        origin_z=float(data.get("origin_z", 0.0)),
        scale=float(data.get("scale", scale)),
    )


def model_point(x: float, y: float, z_m: float, origin: Origin, args: argparse.Namespace) -> list[float]:
    return [
        (x - origin.origin_x) * origin.scale,
        (y - origin.origin_y) * origin.scale,
        (z_m - origin.origin_z) * origin.scale * args.vertical_exaggeration + args.z_offset,
    ]


def normalized_coord(lon: float, lat: float, args: argparse.Namespace) -> tuple[float, float]:
    u = (lon - args.west) / max(args.east - args.west, 1.0e-9)
    v = (lat - args.south) / max(args.north - args.south, 1.0e-9)
    return u, v


def estimated_boundary_z(name: str, lon: float, lat: float, base_z: float, args: argparse.Namespace) -> float:
    if name == "surface":
        return 0.0
    u, v = normalized_coord(lon, lat, args)
    east_deep = u
    bay_deep = max(0.0, 0.75 - v) / 0.75
    upland_shallow = max(0.0, 0.42 - u) * max(0.0, v - 0.35)
    wave = math.sin(u * math.pi * 2.0) * math.cos(v * math.pi * 1.5)
    if name == "base_fill":
        return base_z - 1.0 * bay_deep + 0.3 * wave
    if name == "base_alluvium":
        return base_z - 22.0 * east_deep - 12.0 * bay_deep + 10.0 * upland_shallow + 2.0 * wave
    if name == "base_tokyo_upper":
        return base_z - 28.0 * east_deep - 8.0 * bay_deep + 12.0 * upland_shallow + 3.0 * wave
    if name == "base_tokyo_lower":
        return base_z - 38.0 * east_deep - 10.0 * bay_deep + 15.0 * upland_shallow + 4.0 * wave
    if name == "base_kazusa_top":
        return base_z - 55.0 * east_deep - 18.0 * bay_deep + 22.0 * upland_shallow + 6.0 * wave
    return base_z


def load_surface_points(path: str) -> dict[str, list[tuple[float, float, float]]]:
    if not path:
        return {}
    points: dict[str, list[tuple[float, float, float]]] = {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = row.get("boundary") or row.get("layer") or row.get("name")
            if not name:
                continue
            lon = float(row["lon"])
            lat = float(row["lat"])
            z_m = float(row["z_m"])
            points.setdefault(name, []).append((lon, lat, z_m))
    return points


def idw_z(lon: float, lat: float, points: list[tuple[float, float, float]], fallback: float) -> float:
    if not points:
        return fallback
    total_w = 0.0
    total_z = 0.0
    for p_lon, p_lat, p_z in points:
        d = math.hypot((lon - p_lon) * 90000.0, (lat - p_lat) * 111000.0)
        if d < 1.0:
            return p_z
        w = 1.0 / (d * d)
        total_w += w
        total_z += p_z * w
    return total_z / total_w if total_w > 0.0 else fallback


def make_grid(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lons = np.linspace(args.west, args.east, args.grid_x)
    lats = np.linspace(args.south, args.north, args.grid_y)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    transformer = Transformer.from_crs(WGS84, TOKYO_METERS, always_xy=True)
    x_grid, y_grid = transformer.transform(lon_grid, lat_grid)
    return lon_grid, lat_grid, np.asarray(x_grid), np.asarray(y_grid)


def build_boundaries(args: argparse.Namespace) -> dict[str, np.ndarray]:
    lon_grid, lat_grid, x_grid, y_grid = make_grid(args)
    csv_points = load_surface_points(args.surface_csv)
    z_by_name: dict[str, np.ndarray] = {}
    for name, base_z in DEFAULT_BOUNDARIES:
        z = np.zeros_like(lon_grid, dtype=float)
        points = csv_points.get(name, [])
        for row in range(z.shape[0]):
            for col in range(z.shape[1]):
                fallback = estimated_boundary_z(name, float(lon_grid[row, col]), float(lat_grid[row, col]), base_z, args)
                z[row, col] = idw_z(float(lon_grid[row, col]), float(lat_grid[row, col]), points, fallback)
        z_by_name[name] = z
    z_by_name["_x"] = x_grid
    z_by_name["_y"] = y_grid
    return z_by_name


def vertices_for_surface(x_grid: np.ndarray, y_grid: np.ndarray, z_grid: np.ndarray, origin: Origin, args: argparse.Namespace) -> list[list[float]]:
    rows, cols = z_grid.shape
    vertices = []
    for row in range(rows):
        for col in range(cols):
            vertices.append(model_point(float(x_grid[row, col]), float(y_grid[row, col]), float(z_grid[row, col]), origin, args))
    return vertices


def grid_faces(rows: int, cols: int, reverse: bool = False, offset: int = 0) -> list[list[int]]:
    faces = []
    for row in range(rows - 1):
        for col in range(cols - 1):
            a = offset + row * cols + col
            b = offset + row * cols + col + 1
            c = offset + (row + 1) * cols + col
            d = offset + (row + 1) * cols + col + 1
            if reverse:
                faces.append([a, c, b])
                faces.append([b, c, d])
            else:
                faces.append([a, b, c])
                faces.append([b, d, c])
    return faces


def make_surface_mesh(x_grid: np.ndarray, y_grid: np.ndarray, z_grid: np.ndarray, origin: Origin, args: argparse.Namespace) -> trimesh.Trimesh:
    rows, cols = z_grid.shape
    vertices = vertices_for_surface(x_grid, y_grid, z_grid, origin, args)
    faces = grid_faces(rows, cols)
    return trimesh.Trimesh(vertices=np.asarray(vertices), faces=np.asarray(faces), process=False)


def add_side_faces(faces: list[list[int]], rows: int, cols: int, bottom_offset: int) -> None:
    for col in range(cols - 1):
        a = col
        b = col + 1
        c = bottom_offset + col
        d = bottom_offset + col + 1
        faces.append([a, c, b])
        faces.append([b, c, d])
        a2 = (rows - 1) * cols + col
        b2 = (rows - 1) * cols + col + 1
        c2 = bottom_offset + (rows - 1) * cols + col
        d2 = bottom_offset + (rows - 1) * cols + col + 1
        faces.append([a2, b2, c2])
        faces.append([b2, d2, c2])
    for row in range(rows - 1):
        a = row * cols
        b = (row + 1) * cols
        c = bottom_offset + row * cols
        d = bottom_offset + (row + 1) * cols
        faces.append([a, b, c])
        faces.append([b, d, c])
        a2 = row * cols + cols - 1
        b2 = (row + 1) * cols + cols - 1
        c2 = bottom_offset + row * cols + cols - 1
        d2 = bottom_offset + (row + 1) * cols + cols - 1
        faces.append([a2, c2, b2])
        faces.append([b2, c2, d2])


def make_layer_mesh(
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    top_z: np.ndarray,
    bottom_z: np.ndarray,
    origin: Origin,
    args: argparse.Namespace,
) -> trimesh.Trimesh:
    rows, cols = top_z.shape
    top_vertices = vertices_for_surface(x_grid, y_grid, top_z, origin, args)
    bottom_vertices = vertices_for_surface(x_grid, y_grid, bottom_z, origin, args)
    vertices = top_vertices + bottom_vertices
    bottom_offset = len(top_vertices)
    faces = grid_faces(rows, cols, reverse=False, offset=0)
    faces.extend(grid_faces(rows, cols, reverse=True, offset=bottom_offset))
    add_side_faces(faces, rows, cols, bottom_offset)
    return trimesh.Trimesh(vertices=np.asarray(vertices), faces=np.asarray(faces), process=False)


def export_mesh(mesh: trimesh.Trimesh | None, path: Path) -> None:
    if mesh is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(path)


def write_boundary_csv(path: Path, z_by_name: dict[str, np.ndarray], args: argparse.Namespace) -> None:
    lon_grid, lat_grid, _x_grid, _y_grid = make_grid(args)
    names = [name for name, _base in DEFAULT_BOUNDARIES]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["boundary", "lon", "lat", "z_m"])
        for name in names:
            z = z_by_name[name]
            for row in range(z.shape[0]):
                for col in range(z.shape[1]):
                    writer.writerow([name, lon_grid[row, col], lat_grid[row, col], z[row, col]])


def build(args: argparse.Namespace) -> int:
    origin = load_origin(Path(args.origin_json), args.scale)
    z_by_name = build_boundaries(args)
    x_grid = z_by_name["_x"]
    y_grid = z_by_name["_y"]
    boundary_meshes = []
    layer_meshes = []

    for name, _base_z in DEFAULT_BOUNDARIES:
        mesh = make_surface_mesh(x_grid, y_grid, z_by_name[name], origin, args)
        boundary_meshes.append(mesh)
        export_mesh(mesh, Path(args.out_dir) / f"geology_boundary_{name}.stl")

    for layer_name, top_name, bottom_name in LAYER_PAIRS:
        mesh = make_layer_mesh(x_grid, y_grid, z_by_name[top_name], z_by_name[bottom_name], origin, args)
        layer_meshes.append(mesh)
        export_mesh(mesh, Path(args.out_dir) / f"geology_layer_{layer_name}.stl")

    export_mesh(trimesh.util.concatenate(boundary_meshes), Path(args.out_dir) / args.boundaries_out)
    export_mesh(trimesh.util.concatenate(layer_meshes), Path(args.out_dir) / args.layers_out)
    write_boundary_csv(Path(args.out_dir) / args.sample_csv_out, z_by_name, args)

    print(f"grid: {args.grid_x} x {args.grid_y}")
    print(f"boundaries: {len(boundary_meshes)}")
    print(f"layers: {len(layer_meshes)}")
    print(f"vertical_exaggeration: {args.vertical_exaggeration}")
    print(f"output_dir: {args.out_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build rough geology layer STL surfaces.")
    parser.add_argument("--origin-json", default="origin_epsg6677_zero.json")
    parser.add_argument("--out-dir", default="output/geology_tokyo23_z50")
    parser.add_argument("--surface-csv", default="")
    parser.add_argument("--boundaries-out", default="geology_boundaries_combined.stl")
    parser.add_argument("--layers-out", default="geology_layers_combined.stl")
    parser.add_argument("--sample-csv-out", default="geology_boundary_samples.csv")
    parser.add_argument("--west", type=float, default=139.56)
    parser.add_argument("--east", type=float, default=139.93)
    parser.add_argument("--south", type=float, default=35.52)
    parser.add_argument("--north", type=float, default=35.82)
    parser.add_argument("--grid-x", type=int, default=64)
    parser.add_argument("--grid-y", type=int, default=52)
    parser.add_argument("--scale", type=float, default=0.0005)
    parser.add_argument("--vertical-exaggeration", type=float, default=50.0)
    parser.add_argument("--z-offset", type=float, default=0.0)
    return parser


def main() -> int:
    return build(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

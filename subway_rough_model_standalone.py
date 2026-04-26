#!/usr/bin/env python3
"""Standalone rough subway STL builder with meter-based radii."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TOKYO_METERS = "EPSG:6677"


@dataclass
class Origin:
    origin_x: float
    origin_y: float
    origin_z: float
    scale: float
    crs: str = TOKYO_METERS


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_origin(path: Path | None, rows: list[dict[str, str]], scale: float) -> Origin:
    if path and path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        return Origin(
            origin_x=float(data.get("origin_x", 0.0)),
            origin_y=float(data.get("origin_y", 0.0)),
            origin_z=float(data.get("origin_z", 0.0)),
            scale=float(data.get("scale", scale)),
            crs=str(data.get("crs", TOKYO_METERS)),
        )
    xs = [float(row["x_epsg6677"]) for row in rows if row.get("x_epsg6677")]
    ys = [float(row["y_epsg6677"]) for row in rows if row.get("y_epsg6677")]
    if not xs or not ys:
        return Origin(0.0, 0.0, 0.0, scale)
    return Origin((min(xs) + max(xs)) * 0.5, (min(ys) + max(ys)) * 0.5, 0.0, scale)


def model_point(row: dict[str, str], origin: Origin, z_key: str, vertical_exaggeration: float, z_offset: float) -> list[float]:
    return [
        (float(row["x_epsg6677"]) - origin.origin_x) * origin.scale,
        (float(row["y_epsg6677"]) - origin.origin_y) * origin.scale,
        (float(row[z_key]) - origin.origin_z) * origin.scale * vertical_exaggeration + z_offset,
    ]


def normalize_name(value: str) -> str:
    out = "".join(ch.lower() if ch.isalnum() else "_" for ch in value)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def concatenate_meshes(meshes: list[Any]) -> Any:
    if not meshes:
        return None
    import trimesh

    return trimesh.util.concatenate(meshes)


def export_mesh(mesh: Any, path: Path) -> None:
    if mesh is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(path)


def make_cylinder(p1: list[float], p2: list[float], radius: float, sides: int) -> Any | None:
    import numpy as np
    import trimesh

    if np.linalg.norm(np.asarray(p2) - np.asarray(p1)) <= 0.0:
        return None
    return trimesh.creation.cylinder(radius=radius, sections=sides, segment=[p1, p2])


def build(args: argparse.Namespace) -> int:
    import trimesh  # noqa: F401

    normalized_dir = Path(args.normalized_dir)
    station_rows = read_csv(normalized_dir / "normalized_subway_stations_epsg6677.csv")
    track_rows = read_csv(normalized_dir / "normalized_subway_track_points_epsg6677.csv")
    origin = read_origin(Path(args.origin_json) if args.origin_json else None, station_rows + track_rows, args.scale)

    tube_radius = args.tube_radius_m * origin.scale
    shaft_radius = args.shaft_radius_m * origin.scale
    line_meshes: dict[tuple[str, str], list[Any]] = defaultdict(list)
    all_tube_meshes = []

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in track_rows:
        key = (row["canonical_operator"], row["canonical_line"], row["segment_id"])
        grouped[key].append(row)

    for (operator, line, _segment_id), rows in grouped.items():
        rows.sort(key=lambda row: int(row["sample_index"]))
        for idx in range(len(rows) - 1):
            p1 = model_point(rows[idx], origin, "track_z_m", args.vertical_exaggeration, args.z_offset)
            p2 = model_point(rows[idx + 1], origin, "track_z_m", args.vertical_exaggeration, args.z_offset)
            mesh = make_cylinder(p1, p2, tube_radius, args.tube_sides)
            if mesh is None:
                continue
            all_tube_meshes.append(mesh)
            line_meshes[(operator, line)].append(mesh)

    shaft_meshes = []
    if args.station_shafts:
        for row in station_rows:
            p1 = model_point(row, origin, "station_z_m", args.vertical_exaggeration, args.z_offset)
            surface = dict(row)
            surface["surface_z_m"] = row.get("surface_z_m") or "0"
            p2 = model_point(surface, origin, "surface_z_m", args.vertical_exaggeration, args.z_offset)
            mesh = make_cylinder(p1, p2, shaft_radius, args.shaft_sides)
            if mesh is not None:
                shaft_meshes.append(mesh)

    tube_mesh = concatenate_meshes(all_tube_meshes)
    shaft_mesh = concatenate_meshes(shaft_meshes)
    combined = concatenate_meshes([mesh for mesh in [tube_mesh, shaft_mesh] if mesh is not None])

    export_mesh(tube_mesh, Path(args.tube_out))
    export_mesh(shaft_mesh, Path(args.shaft_out))
    export_mesh(combined, Path(args.out))

    if args.split_by_line:
        by_line_dir = Path(args.by_line_dir)
        by_line_dir.mkdir(parents=True, exist_ok=True)
        for (operator, line), meshes in line_meshes.items():
            name = normalize_name(f"{operator}_{line}")
            export_mesh(concatenate_meshes(meshes), by_line_dir / f"{name}.stl")

    print(f"origin_x: {origin.origin_x}")
    print(f"origin_y: {origin.origin_y}")
    print(f"scale: {origin.scale}")
    print(f"tube_radius_model: {tube_radius}")
    print(f"shaft_radius_model: {shaft_radius}")
    print(f"tube segments: {len(all_tube_meshes)}")
    print(f"shafts: {len(shaft_meshes)}")
    print(f"output: {args.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build standalone rough subway STL with meter-based radii.")
    parser.add_argument("--normalized-dir", required=True)
    parser.add_argument("--origin-json", default=None)
    parser.add_argument("--out", default="output/subway_rough_with_shafts_mscale.stl")
    parser.add_argument("--tube-out", default="output/subway_rough_all_lines_mscale.stl")
    parser.add_argument("--shaft-out", default="output/subway_station_shafts_mscale.stl")
    parser.add_argument("--by-line-dir", default="output/subway_rough_by_line_mscale")
    parser.add_argument("--scale", type=float, default=0.0005)
    parser.add_argument("--tube-radius-m", type=float, default=10.0)
    parser.add_argument("--shaft-radius-m", type=float, default=16.0)
    parser.add_argument("--tube-sides", type=int, default=8)
    parser.add_argument("--shaft-sides", type=int, default=8)
    parser.add_argument("--vertical-exaggeration", type=float, default=2.0)
    parser.add_argument("--z-offset", type=float, default=0.0)
    parser.add_argument("--station-shafts", action="store_true")
    parser.add_argument("--split-by-line", action="store_true")
    return parser


def main() -> int:
    return build(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

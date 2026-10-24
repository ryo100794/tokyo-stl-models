"""Subway normalization and rough STL generation helpers."""

from __future__ import annotations

import csv
import json
import math
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from pyproj import Transformer


WGS84 = "EPSG:4326"
TOKYO_METERS = "EPSG:6677"

DEFAULT_BBOX = {
    "west": 139.655,
    "east": 139.735,
    "south": 35.640,
    "north": 35.715,
}

STATION_COLUMNS = [
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

TRACK_POINT_COLUMNS = [
    "canonical_operator",
    "canonical_line",
    "segment_id",
    "sample_index",
    "lon",
    "lat",
    "x_epsg6677",
    "y_epsg6677",
    "surface_z_m",
    "depth_m",
    "track_z_m",
    "interpolation_method",
    "confidence",
    "note",
]

DEPTH_COLUMNS = [
    "operator",
    "line",
    "station",
    "station_id",
    "lon",
    "lat",
    "depth_m",
    "depth_type",
    "source_url",
    "source_label",
    "confidence",
    "note",
]

ALIAS_COLUMNS = [
    "canonical_operator",
    "canonical_line",
    "canonical_station",
    "alias_operator",
    "alias_line",
    "alias_station",
    "source",
]

CONFIDENCE_RANK = {
    "official": 0,
    "literature": 1,
    "secondary": 2,
    "estimated": 3,
    "unknown": 4,
}

ROUTES = [
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Ginza Line",
        "default_depth_m": 12.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Ginza Line", "Tokyo Metro Ginza Line", "Ginza", "G"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Marunouchi Line",
        "default_depth_m": 18.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Marunouchi Line", "Tokyo Metro Marunouchi Line", "Marunouchi", "M"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Hibiya Line",
        "default_depth_m": 18.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Hibiya Line", "Tokyo Metro Hibiya Line", "Hibiya", "H"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Tozai Line",
        "default_depth_m": 20.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Tozai Line", "Tokyo Metro Tozai Line", "Tozai", "T"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Chiyoda Line",
        "default_depth_m": 23.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Chiyoda Line", "Tokyo Metro Chiyoda Line", "Chiyoda", "C"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Yurakucho Line",
        "default_depth_m": 25.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Yurakucho Line", "Tokyo Metro Yurakucho Line", "Yurakucho", "Y"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Hanzomon Line",
        "default_depth_m": 28.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Hanzomon Line", "Tokyo Metro Hanzomon Line", "Hanzomon", "Z"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Namboku Line",
        "default_depth_m": 28.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Namboku Line", "Tokyo Metro Namboku Line", "Namboku", "N"],
        },
    },
    {
        "canonical_operator": "Tokyo Metro",
        "canonical_line": "Fukutoshin Line",
        "default_depth_m": 30.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Tokyo Metro", "Tokyo Subway", "Tokyo Metro Co Ltd"],
            "line": ["Fukutoshin Line", "Tokyo Metro Fukutoshin Line", "Fukutoshin", "F"],
        },
    },
    {
        "canonical_operator": "Toei Subway",
        "canonical_line": "Asakusa Line",
        "default_depth_m": 20.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Toei Subway", "Tokyo Metropolitan Bureau of Transportation", "Toei"],
            "line": ["Asakusa Line", "Toei Asakusa Line", "Asakusa", "A"],
        },
    },
    {
        "canonical_operator": "Toei Subway",
        "canonical_line": "Mita Line",
        "default_depth_m": 24.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Toei Subway", "Tokyo Metropolitan Bureau of Transportation", "Toei"],
            "line": ["Mita Line", "Toei Mita Line", "Mita", "I"],
        },
    },
    {
        "canonical_operator": "Toei Subway",
        "canonical_line": "Shinjuku Line",
        "default_depth_m": 25.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Toei Subway", "Tokyo Metropolitan Bureau of Transportation", "Toei"],
            "line": ["Shinjuku Line", "Toei Shinjuku Line", "Shinjuku", "S"],
        },
    },
    {
        "canonical_operator": "Toei Subway",
        "canonical_line": "Oedo Line",
        "default_depth_m": 35.0,
        "tube_radius": 1.0,
        "aliases": {
            "operator": ["Toei Subway", "Tokyo Metropolitan Bureau of Transportation", "Toei"],
            "line": ["Oedo Line", "Toei Oedo Line", "Oedo", "E"],
        },
    },
]

SAMPLE_STATIONS = [
    ("Tokyo Metro", "Ginza Line", "Shibuya", "sample_ginza_shibuya", 139.7012, 35.6580, 12.0),
    ("Tokyo Metro", "Ginza Line", "Omote-sando", "sample_ginza_omotesando", 139.7122, 35.6652, 14.0),
    ("Tokyo Metro", "Ginza Line", "Aoyama-itchome", "sample_ginza_aoyama", 139.7240, 35.6728, 15.0),
    ("Tokyo Metro", "Marunouchi Line", "Shinjuku", "sample_m_shinjuku", 139.7005, 35.6906, 18.0),
    ("Tokyo Metro", "Marunouchi Line", "Shinjuku-sanchome", "sample_m_sanchome", 139.7065, 35.6901, 20.0),
    ("Tokyo Metro", "Marunouchi Line", "Yotsuya-sanchome", "sample_m_yotsuya3", 139.7205, 35.6879, 18.0),
    ("Tokyo Metro", "Hibiya Line", "Ebisu", "sample_h_ebisu", 139.7101, 35.6467, 18.0),
    ("Tokyo Metro", "Hibiya Line", "Hiro-o", "sample_h_hiroo", 139.7223, 35.6516, 20.0),
    ("Tokyo Metro", "Tozai Line", "Kagurazaka", "sample_t_kagurazaka", 139.7350, 35.7038, 20.0),
    ("Tokyo Metro", "Chiyoda Line", "Yoyogi-uehara", "sample_c_yoyogiuehara", 139.6804, 35.6690, 20.0),
    ("Tokyo Metro", "Chiyoda Line", "Yoyogi-koen", "sample_c_yoyogikoen", 139.6891, 35.6690, 22.0),
    ("Tokyo Metro", "Chiyoda Line", "Meiji-jingumae", "sample_c_meiji", 139.7055, 35.6691, 24.0),
    ("Tokyo Metro", "Yurakucho Line", "Iidabashi", "sample_y_iidabashi", 139.7448, 35.7021, 25.0),
    ("Tokyo Metro", "Hanzomon Line", "Shibuya", "sample_z_shibuya", 139.7012, 35.6580, 28.0),
    ("Tokyo Metro", "Hanzomon Line", "Omote-sando", "sample_z_omotesando", 139.7122, 35.6652, 28.0),
    ("Tokyo Metro", "Hanzomon Line", "Aoyama-itchome", "sample_z_aoyama", 139.7240, 35.6728, 30.0),
    ("Tokyo Metro", "Namboku Line", "Yotsuya", "sample_n_yotsuya", 139.7300, 35.6860, 28.0),
    ("Tokyo Metro", "Fukutoshin Line", "Shibuya", "sample_f_shibuya", 139.7012, 35.6580, 30.0),
    ("Tokyo Metro", "Fukutoshin Line", "Meiji-jingumae", "sample_f_meiji", 139.7055, 35.6691, 30.0),
    ("Tokyo Metro", "Fukutoshin Line", "Kita-sando", "sample_f_kitasando", 139.7056, 35.6787, 30.0),
    ("Tokyo Metro", "Fukutoshin Line", "Shinjuku-sanchome", "sample_f_sanchome", 139.7065, 35.6901, 30.0),
    ("Toei Subway", "Asakusa Line", "Daimon", "sample_a_daimon", 139.7546, 35.6566, 20.0),
    ("Toei Subway", "Mita Line", "Onarimon", "sample_i_onarimon", 139.7516, 35.6603, 24.0),
    ("Toei Subway", "Shinjuku Line", "Shinjuku", "sample_s_shinjuku", 139.7005, 35.6906, 25.0),
    ("Toei Subway", "Shinjuku Line", "Shinjuku-sanchome", "sample_s_sanchome", 139.7065, 35.6901, 25.0),
    ("Toei Subway", "Shinjuku Line", "Akebonobashi", "sample_s_akebonobashi", 139.7228, 35.6924, 27.0),
    ("Toei Subway", "Oedo Line", "Yoyogi", "sample_e_yoyogi", 139.7023, 35.6830, 35.0),
    ("Toei Subway", "Oedo Line", "Shinjuku", "sample_e_shinjuku", 139.7005, 35.6906, 38.0),
    ("Toei Subway", "Oedo Line", "Tochomae", "sample_e_tochomae", 139.6922, 35.6906, 38.0),
]


@dataclass
class Origin:
    origin_x: float
    origin_y: float
    origin_z: float
    scale: float
    crs: str = TOKYO_METERS


def read_yaml(path: Path) -> Any:
    import yaml

    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def write_yaml(path: Path, data: Any) -> None:
    import yaml

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        yaml.safe_dump(data, fh, sort_keys=False, allow_unicode=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def command_subway_init(args: Any) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_yaml(out_dir / "subway_lines.yaml", {"routes": ROUTES})
    write_yaml(
        out_dir / "source_priority.yaml",
        {
            "position_priority": ["manual_csv", "gtfs", "n02", "projected_to_line"],
            "depth_priority": ["official", "literature", "secondary", "estimated", "unknown"],
        },
    )
    depth_rows = []
    if args.sample:
        for op, line, station, station_id, lon, lat, depth in SAMPLE_STATIONS:
            depth_rows.append(
                {
                    "operator": op,
                    "line": line,
                    "station": station,
                    "station_id": station_id,
                    "lon": lon,
                    "lat": lat,
                    "depth_m": depth,
                    "depth_type": "estimated",
                    "source_url": "",
                    "source_label": "sample_seed",
                    "confidence": "estimated",
                    "note": "Approximate seed for visualization; replace with verified source.",
                }
            )
    write_csv(out_dir / "subway_depth.csv", depth_rows, DEPTH_COLUMNS)
    write_csv(out_dir / "station_aliases.csv", [], ALIAS_COLUMNS)
    (out_dir / "README_subway_data.md").write_text(
        "\n".join(
            [
                "# Subway Data Inputs",
                "",
                "Use public railway data, GTFS/GTFS-JP, and verified manual depth CSV rows.",
                "The sample rows are approximate visualization seeds.",
                "Replace sample depths with verified values before analytical use.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"created: {out_dir}")
    return 0


def bbox_contains(row: dict[str, Any], bbox: dict[str, float]) -> bool:
    lon = float(row["lon"])
    lat = float(row["lat"])
    return bbox["west"] <= lon <= bbox["east"] and bbox["south"] <= lat <= bbox["north"]


def make_alias_index(routes: list[dict[str, Any]]) -> dict[tuple[str, str], tuple[str, str]]:
    out: dict[tuple[str, str], tuple[str, str]] = {}
    for route in routes:
        cop = route["canonical_operator"]
        cline = route["canonical_line"]
        op_values = [cop] + route.get("aliases", {}).get("operator", [])
        line_values = [cline] + route.get("aliases", {}).get("line", [])
        for op in op_values:
            for line in line_values:
                out[(normalize_key(op), normalize_key(line))] = (cop, cline)
    return out


def normalize_key(value: Any) -> str:
    return str(value or "").strip().lower().replace("_", " ").replace("-", " ")


def route_defaults(routes: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, Any]]:
    return {(r["canonical_operator"], r["canonical_line"]): r for r in routes}


def canonicalize_row(row: dict[str, str], alias_index: dict[tuple[str, str], tuple[str, str]]) -> tuple[str, str] | None:
    key = (normalize_key(row.get("operator")), normalize_key(row.get("line")))
    return alias_index.get(key)


def safe_float(value: Any, default: float | None = None) -> float | None:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except ValueError:
        return default


def load_manual_stations(depth_csv: Path, routes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    alias_index = make_alias_index(routes)
    defaults = route_defaults(routes)
    rows: list[dict[str, Any]] = []
    for row in read_csv(depth_csv):
        canon = canonicalize_row(row, alias_index)
        if canon is None:
            continue
        cop, cline = canon
        default_depth = float(defaults[(cop, cline)].get("default_depth_m", 25.0))
        depth = safe_float(row.get("depth_m"), default_depth)
        lon = safe_float(row.get("lon"))
        lat = safe_float(row.get("lat"))
        if lon is None or lat is None:
            continue
        confidence = row.get("confidence") or "unknown"
        depth_type = row.get("depth_type") or "unknown"
        note = row.get("note") or ""
        if row.get("depth_m") in ("", None):
            confidence = "estimated"
            depth_type = "estimated"
            note = append_note(note, "default route depth used")
        rows.append(
            {
                "canonical_operator": cop,
                "canonical_line": cline,
                "canonical_station": row.get("station") or row.get("station_id") or "unknown",
                "station_id": row.get("station_id") or make_station_id(cop, cline, row.get("station")),
                "lon": lon,
                "lat": lat,
                "surface_z_m": 0.0,
                "depth_m": depth,
                "station_z_m": -float(depth or 0.0),
                "depth_type": depth_type,
                "confidence": confidence,
                "source_label": row.get("source_label") or "manual_csv",
                "source_url": row.get("source_url") or "",
                "note": note,
                "source_candidates": [{"source": "manual_csv", "row": row}],
            }
        )
    return rows


def read_gtfs_table(gtfs_path: Path, name: str) -> list[dict[str, str]]:
    if gtfs_path.is_dir():
        return read_csv(gtfs_path / name)
    if gtfs_path.is_file() and gtfs_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(gtfs_path, "r") as zf:
            try:
                with zf.open(name, "r") as raw:
                    text = raw.read().decode("utf-8-sig")
            except KeyError:
                return []
        return list(csv.DictReader(text.splitlines()))
    return []


def load_gtfs_stations(gtfs_dirs: list[str], routes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    alias_index = make_alias_index(routes)
    defaults = route_defaults(routes)
    out: list[dict[str, Any]] = []
    for item in gtfs_dirs:
        gtfs_path = Path(item)
        if not gtfs_path.exists():
            continue
        gtfs_routes = read_gtfs_table(gtfs_path, "routes.txt")
        stops = {row.get("stop_id", ""): row for row in read_gtfs_table(gtfs_path, "stops.txt")}
        trips = read_gtfs_table(gtfs_path, "trips.txt")
        stop_times = read_gtfs_table(gtfs_path, "stop_times.txt")
        route_lookup: dict[str, tuple[str, str]] = {}
        for row in gtfs_routes:
            name_values = [row.get("route_long_name"), row.get("route_short_name")]
            for name in name_values:
                for (alias_op, alias_line), canon in alias_index.items():
                    if normalize_key(name) == alias_line:
                        route_lookup[row.get("route_id", "")] = canon
                        break
                if row.get("route_id", "") in route_lookup:
                    break
        trip_to_route = {row.get("trip_id", ""): row.get("route_id", "") for row in trips}
        route_stop_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in stop_times:
            route_id = trip_to_route.get(row.get("trip_id", ""))
            canon = route_lookup.get(route_id)
            if canon is None:
                continue
            stop = stops.get(row.get("stop_id", ""))
            if not stop:
                continue
            merged = dict(stop)
            merged["stop_sequence"] = row.get("stop_sequence", "0")
            route_stop_rows[canon].append(merged)
        for (cop, cline), rows in route_stop_rows.items():
            seen = set()
            rows.sort(key=lambda row: int(row.get("stop_sequence") or 0))
            default_depth = float(defaults[(cop, cline)].get("default_depth_m", 25.0))
            for stop in rows:
                stop_id = stop.get("stop_id", "")
                if stop_id in seen:
                    continue
                seen.add(stop_id)
                lon = safe_float(stop.get("stop_lon"))
                lat = safe_float(stop.get("stop_lat"))
                if lon is None or lat is None:
                    continue
                name = stop.get("stop_name") or stop_id
                out.append(
                    {
                        "canonical_operator": cop,
                        "canonical_line": cline,
                        "canonical_station": name,
                        "station_id": stop_id,
                        "lon": lon,
                        "lat": lat,
                        "surface_z_m": 0.0,
                        "depth_m": default_depth,
                        "station_z_m": -default_depth,
                        "depth_type": "estimated",
                        "confidence": "estimated",
                        "source_label": "gtfs",
                        "source_url": "",
                        "note": "GTFS position with route default depth.",
                        "source_candidates": [{"source": "gtfs", "row": stop}],
                    }
                )
    return out


def merge_station_candidates(primary: list[dict[str, Any]], secondary: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = list(primary)
    keys = {station_identity(row) for row in out}
    for row in secondary:
        key = station_identity(row)
        if key in keys:
            continue
        keys.add(key)
        out.append(row)
    return out


def station_identity(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        normalize_key(row.get("canonical_operator")),
        normalize_key(row.get("canonical_line")),
        normalize_key(row.get("canonical_station")),
    )


def append_note(note: str, extra: str) -> str:
    if not note:
        return extra
    return f"{note}; {extra}"


def make_station_id(operator: str, line: str, station: Any) -> str:
    base = "_".join([operator, line, str(station or "unknown")])
    return normalize_key(base).replace(" ", "_")


def add_epsg6677(rows: list[dict[str, Any]]) -> None:
    transformer = Transformer.from_crs(WGS84, TOKYO_METERS, always_xy=True)
    for row in rows:
        x, y = transformer.transform(float(row["lon"]), float(row["lat"]))
        row["x_epsg6677"] = x
        row["y_epsg6677"] = y


def filter_rows(rows: list[dict[str, Any]], bbox: dict[str, float] | None, clip_mode: str) -> list[dict[str, Any]]:
    if clip_mode == "all" or bbox is None:
        return rows
    return [row for row in rows if bbox_contains(row, bbox)]


def station_sort_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (row["canonical_operator"], row["canonical_line"], str(row.get("station_id") or row.get("canonical_station")))


def build_lines(stations: list[dict[str, Any]], clip_mode: str) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in sorted(stations, key=station_sort_key):
        grouped[(row["canonical_operator"], row["canonical_line"])].append(row)
    lines = []
    for (op, line), rows in grouped.items():
        if len(rows) < 2:
            continue
        coords = [[float(row["lon"]), float(row["lat"])] for row in rows]
        lines.append(
            {
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": coords},
                "properties": {
                    "canonical_operator": op,
                    "canonical_line": line,
                    "line_id": make_station_id(op, line, "line"),
                    "source_priority": "manual_csv",
                    "clip_mode": clip_mode,
                    "confidence": merged_confidence(rows),
                    "note": "Rough line from station sequence.",
                },
            }
        )
    return lines


def load_n02_lines(railway_root: str | None, routes: list[dict[str, Any]], bbox: dict[str, float] | None, clip_mode: str) -> list[dict[str, Any]]:
    if not railway_root:
        return []
    root = Path(railway_root)
    if not root.exists():
        return []
    alias_index = make_alias_index(routes)
    features = []
    for path in list(root.rglob("*.geojson")) + list(root.rglob("*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        for feature in data.get("features", []):
            canon = match_feature_route(feature.get("properties", {}), alias_index)
            if canon is None:
                continue
            geom = feature.get("geometry") or {}
            line_parts = geometry_to_lines(geom)
            for coords in line_parts:
                if bbox and clip_mode == "bbox" and not any(coord_in_bbox(coord, bbox) for coord in coords):
                    continue
                features.append(
                    {
                        "type": "Feature",
                        "geometry": {"type": "LineString", "coordinates": coords},
                        "properties": {
                            "canonical_operator": canon[0],
                            "canonical_line": canon[1],
                            "line_id": make_station_id(canon[0], canon[1], path.stem),
                            "source_priority": "n02",
                            "clip_mode": clip_mode,
                            "confidence": "secondary",
                            "note": f"Line geometry from {path.name}.",
                        },
                    }
                )
    return features


def match_feature_route(props: dict[str, Any], alias_index: dict[tuple[str, str], tuple[str, str]]) -> tuple[str, str] | None:
    """Match whole aliases in designated route fields, never arbitrary substrings.

    A known operator disambiguates short route codes. Unknown/conflicting
    operator fields and ambiguous line aliases deliberately remain unresolved.
    N02 Japanese aliases must be supplied in the route configuration just like
    other provider aliases; an unrelated property is never guessed to be a line.
    """
    operator_keys = {"operator", "operator_name", "canonical_operator", "agency_name", "n02_004"}
    line_keys = {"line", "line_name", "canonical_line", "route_name", "route_long_name", "route_short_name", "n02_003"}
    operators = {normalize_key(value) for key, value in props.items()
                 if str(key).lower() in operator_keys and normalize_key(value)}
    lines = {normalize_key(value) for key, value in props.items()
             if str(key).lower() in line_keys and normalize_key(value)}
    candidates = set()
    for line in lines:
        line_matches = {canon for (operator_alias, line_alias), canon in alias_index.items()
                        if line_alias == line and (not operators or operator_alias in operators)}
        candidates.update(line_matches)
    if operators:
        recognized_operators = {operator_alias for operator_alias, _ in alias_index}
        if not operators.issubset(recognized_operators):
            return None
        canonical_operators = {canon[0] for (operator_alias, _), canon in alias_index.items()
                               if operator_alias in operators}
        if len(canonical_operators) != 1:
            return None
    return next(iter(candidates)) if len(candidates) == 1 else None


def geometry_to_lines(geom: dict[str, Any]) -> list[list[list[float]]]:
    gtype = geom.get("type")
    coords = geom.get("coordinates") or []
    if gtype == "LineString":
        return [coords]
    if gtype == "MultiLineString":
        return coords
    return []


def coord_in_bbox(coord: list[float], bbox: dict[str, float]) -> bool:
    lon = float(coord[0])
    lat = float(coord[1])
    return bbox["west"] <= lon <= bbox["east"] and bbox["south"] <= lat <= bbox["north"]


def merged_confidence(rows: list[dict[str, Any]]) -> str:
    """Combined evidence can be no stronger than its least certain input."""
    if not rows:
        return "unknown"
    confidences = [normalize_key(row.get("confidence")) for row in rows]
    confidences = [conf if conf in CONFIDENCE_RANK else "unknown" for conf in confidences]
    return max(confidences, key=CONFIDENCE_RANK.__getitem__)


def build_segments(stations: list[dict[str, Any]], sample_step_m: float = 250.0) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in sorted(stations, key=station_sort_key):
        grouped[(row["canonical_operator"], row["canonical_line"])].append(row)

    features = []
    track_points = []
    for (op, line), rows in grouped.items():
        if len(rows) < 2:
            continue
        for idx in range(len(rows) - 1):
            a = rows[idx]
            b = rows[idx + 1]
            segment_id = f"{make_station_id(op, line, idx)}_{idx + 1}"
            features.append(
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[float(a["lon"]), float(a["lat"])], [float(b["lon"]), float(b["lat"])]],
                    },
                    "properties": {
                        "canonical_operator": op,
                        "canonical_line": line,
                        "segment_id": segment_id,
                        "from_station": a["canonical_station"],
                        "to_station": b["canonical_station"],
                        "interpolation_method": "linear_station_to_station",
                        "confidence": merged_confidence([a, b]),
                        "note": "Rough station-to-station segment.",
                    },
                }
            )
            track_points.extend(sample_segment_points(a, b, segment_id, sample_step_m))
    return features, track_points


def sample_segment_points(a: dict[str, Any], b: dict[str, Any], segment_id: str, sample_step_m: float) -> list[dict[str, Any]]:
    ax = float(a["x_epsg6677"])
    ay = float(a["y_epsg6677"])
    bx = float(b["x_epsg6677"])
    by = float(b["y_epsg6677"])
    dist = math.hypot(bx - ax, by - ay)
    count = max(2, int(math.ceil(dist / sample_step_m)) + 1)
    out = []
    for idx in range(count):
        t = idx / (count - 1)
        row = {
            "canonical_operator": a["canonical_operator"],
            "canonical_line": a["canonical_line"],
            "segment_id": segment_id,
            "sample_index": idx,
            "lon": lerp(float(a["lon"]), float(b["lon"]), t),
            "lat": lerp(float(a["lat"]), float(b["lat"]), t),
            "x_epsg6677": lerp(ax, bx, t),
            "y_epsg6677": lerp(ay, by, t),
            "surface_z_m": lerp(float(a["surface_z_m"]), float(b["surface_z_m"]), t),
            "depth_m": lerp(float(a["depth_m"]), float(b["depth_m"]), t),
            "track_z_m": lerp(float(a["station_z_m"]), float(b["station_z_m"]), t),
            "interpolation_method": "linear_station_to_station",
            "confidence": merged_confidence([a, b]),
            "note": "Sampled from normalized station pair.",
        }
        out.append(row)
    return out


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def write_geojson(path: Path, features: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "type": "FeatureCollection",
        "name": path.stem,
        "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
        "features": features,
    }
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def station_features(stations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    features = []
    for row in stations:
        props = {col: row.get(col, "") for col in STATION_COLUMNS}
        features.append(
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [float(row["lon"]), float(row["lat"])]},
                "properties": props,
            }
        )
    return features


def command_subway_normalize(args: Any) -> int:
    lines_yaml = Path(args.lines_yaml)
    routes = read_yaml(lines_yaml).get("routes", ROUTES)
    manual_stations = load_manual_stations(Path(args.depth_csv), routes)
    gtfs_stations = load_gtfs_stations(args.gtfs_dir, routes)
    stations = merge_station_candidates(manual_stations, gtfs_stations)
    add_epsg6677(stations)
    bbox = make_bbox(args)
    stations = filter_rows(stations, bbox, args.clip_mode)
    stations.sort(key=station_sort_key)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = load_n02_lines(args.railway_root, routes, bbox, args.clip_mode) + build_lines(stations, args.clip_mode)
    segments, track_points = build_segments(stations)

    write_geojson(out_dir / "normalized_subway_stations.geojson", station_features(stations))
    write_geojson(out_dir / "normalized_subway_lines.geojson", lines)
    write_geojson(out_dir / "normalized_subway_segments.geojson", segments)
    write_csv(out_dir / "normalized_subway_stations_epsg6677.csv", stations, STATION_COLUMNS)
    write_csv(out_dir / "normalized_subway_track_points_epsg6677.csv", track_points, TRACK_POINT_COLUMNS)
    write_source_candidates(out_dir / "normalized_subway_source_candidates.json", stations)
    write_report(out_dir, out_dir / "subway_data_quality_report.md")
    print(f"stations: {len(stations)}")
    print(f"lines: {len(lines)}")
    print(f"segments: {len(segments)}")
    print(f"track points: {len(track_points)}")
    print(f"output: {out_dir}")
    return 0


def write_source_candidates(path: Path, stations: list[dict[str, Any]]) -> None:
    rows = []
    for station in stations:
        rows.append(
            {
                "canonical_operator": station.get("canonical_operator"),
                "canonical_line": station.get("canonical_line"),
                "canonical_station": station.get("canonical_station"),
                "station_id": station.get("station_id"),
                "source_candidates": station.get("source_candidates", []),
            }
        )
    path.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")


def make_bbox(args: Any) -> dict[str, float] | None:
    if args.clip_mode == "all":
        return None
    return {"west": args.west, "east": args.east, "south": args.south, "north": args.north}


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


def model_point(row: dict[str, Any], origin: Origin, z_key: str, vertical_exaggeration: float, z_offset: float) -> list[float]:
    return [
        (float(row["x_epsg6677"]) - origin.origin_x) * origin.scale,
        (float(row["y_epsg6677"]) - origin.origin_y) * origin.scale,
        (float(row[z_key]) - origin.origin_z) * origin.scale * vertical_exaggeration + z_offset,
    ]


def command_subway_rough_stl(args: Any) -> int:
    import numpy as np
    import trimesh

    normalized_dir = Path(args.normalized_dir)
    station_rows = read_csv(normalized_dir / "normalized_subway_stations_epsg6677.csv")
    track_rows = read_csv(normalized_dir / "normalized_subway_track_points_epsg6677.csv")
    origin = read_origin(Path(args.origin_json) if args.origin_json else None, station_rows + track_rows, args.scale)
    meshes = []
    line_meshes: dict[tuple[str, str], list[Any]] = defaultdict(list)

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in track_rows:
        grouped[(row["canonical_operator"], row["canonical_line"], row["segment_id"])].append(row)

    for (op, line, _segment_id), rows in grouped.items():
        rows.sort(key=lambda row: int(row["sample_index"]))
        for idx in range(len(rows) - 1):
            p1 = model_point(rows[idx], origin, "track_z_m", args.vertical_exaggeration, args.z_offset)
            p2 = model_point(rows[idx + 1], origin, "track_z_m", args.vertical_exaggeration, args.z_offset)
            if np.linalg.norm(np.asarray(p2) - np.asarray(p1)) <= 0.0:
                continue
            mesh = trimesh.creation.cylinder(radius=args.tube_radius, sections=args.tube_sides, segment=[p1, p2])
            meshes.append(mesh)
            line_meshes[(op, line)].append(mesh)

    tube_mesh = concatenate_meshes(meshes)
    if args.tube_out:
        export_mesh(tube_mesh, Path(args.tube_out))

    shaft_meshes = []
    if args.station_shafts:
        for row in station_rows:
            p1 = model_point(row, origin, "station_z_m", args.vertical_exaggeration, args.z_offset)
            surface_row = dict(row)
            surface_row["surface_z_m"] = row.get("surface_z_m") or "0"
            p2 = model_point(surface_row, origin, "surface_z_m", args.vertical_exaggeration, args.z_offset)
            if np.linalg.norm(np.asarray(p2) - np.asarray(p1)) <= 0.0:
                continue
            shaft_meshes.append(trimesh.creation.cylinder(radius=args.shaft_radius, sections=args.shaft_sides, segment=[p1, p2]))
    shaft_mesh = concatenate_meshes(shaft_meshes)
    if args.shaft_out:
        export_mesh(shaft_mesh, Path(args.shaft_out))

    if args.split_by_line:
        out_dir = Path(args.by_line_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for (op, line), part_meshes in line_meshes.items():
            name = normalize_key(f"{op}_{line}").replace(" ", "_")
            export_mesh(concatenate_meshes(part_meshes), out_dir / f"{name}.stl")

    combined = concatenate_meshes([mesh for mesh in [tube_mesh, shaft_mesh] if mesh is not None])
    export_mesh(combined, Path(args.out))
    if args.debug_points:
        write_debug_points(Path(args.debug_points), station_rows, track_rows, origin, args.vertical_exaggeration, args.z_offset)
    print(f"tube segments: {len(meshes)}")
    print(f"shafts: {len(shaft_meshes)}")
    print(f"output: {args.out}")
    return 0


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


def write_debug_points(
    path: Path,
    station_rows: list[dict[str, str]],
    track_rows: list[dict[str, str]],
    origin: Origin,
    vertical_exaggeration: float,
    z_offset: float,
) -> None:
    points = []
    for row in station_rows:
        points.append(model_point(row, origin, "station_z_m", vertical_exaggeration, z_offset))
    for row in track_rows:
        points.append(model_point(row, origin, "track_z_m", vertical_exaggeration, z_offset))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write("ply\nformat ascii 1.0\n")
        fh.write(f"element vertex {len(points)}\n")
        fh.write("property float x\nproperty float y\nproperty float z\nend_header\n")
        for point in points:
            fh.write(f"{point[0]} {point[1]} {point[2]}\n")


def write_report(normalized_dir: Path, out_path: Path) -> None:
    stations = read_csv(normalized_dir / "normalized_subway_stations_epsg6677.csv")
    track_points = read_csv(normalized_dir / "normalized_subway_track_points_epsg6677.csv")
    by_line: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in stations:
        by_line[(row["canonical_operator"], row["canonical_line"])].append(row)
    lines = ["# Subway Data Quality Report", ""]
    lines.append(f"Station count: {len(stations)}")
    lines.append(f"Track sample point count: {len(track_points)}")
    lines.append("")
    for key in sorted(by_line):
        rows = by_line[key]
        conf_counts: dict[str, int] = defaultdict(int)
        depth_count = 0
        for row in rows:
            conf_counts[row.get("confidence") or "unknown"] += 1
            if row.get("depth_m") not in ("", None):
                depth_count += 1
        lines.append(f"## {key[0]} {key[1]}")
        lines.append("")
        lines.append(f"- stations: {len(rows)}")
        lines.append(f"- stations with depth: {depth_count}")
        for conf in ["official", "literature", "secondary", "estimated", "unknown"]:
            lines.append(f"- {conf}: {conf_counts.get(conf, 0)}")
        estimated = [row["canonical_station"] for row in rows if row.get("confidence") == "estimated"]
        if estimated:
            lines.append(f"- estimated stations: {', '.join(estimated)}")
        lines.append("- adopted source: manual_csv/seed unless replaced by GTFS or N02 in a later run")
        lines.append("")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def command_subway_report(args: Any) -> int:
    write_report(Path(args.normalized_dir), Path(args.out))
    print(f"report: {args.out}")
    return 0

# Standalone N02 Subway Rough Model

This note describes `subway_n02_rough_model_standalone.py`.

The file is independent from the main CLI to avoid conflicts with parallel edits.

## Purpose

The earlier sample-based subway model connected stations with straight lines and had too few stations.

This standalone script uses National Land Numerical Information N02 railway data:

- `UTF-8/N02-24_RailroadSection.geojson` for curved railway line geometry
- `UTF-8/N02-24_Station.geojson` for station locations

It builds rough subway tubes and station shafts for the 13 Tokyo subway lines.

## Tokyo 23 wards command

The current Tokyo 23 wards rough model was generated with:

```bash
/root/.venvs/tokyo-stl-models/bin/python subway_n02_rough_model_standalone.py \
  --n02-zip data/n02/N02-24_GML.zip \
  --normalized-out-dir normalized_subway_n02_tokyo23 \
  --out output/subway_n02_tokyo23_rough_with_shafts_round.stl \
  --tube-out output/subway_n02_tokyo23_rough_all_lines_round.stl \
  --shaft-out output/subway_n02_tokyo23_station_shafts_round.stl \
  --by-line-dir output/subway_n02_tokyo23_by_line_round \
  --west 139.56 \
  --east 139.93 \
  --south 35.52 \
  --north 35.82 \
  --clip-mode bbox \
  --scale 0.0005 \
  --tube-radius-m 10.0 \
  --shaft-radius-m 16.0 \
  --tube-sides 16 \
  --shaft-sides 16 \
  --min-point-spacing-m 50.0 \
  --vertical-exaggeration 2.0 \
  --station-shafts \
  --split-by-line
```

## Current local outputs

These files are ignored by Git:

```text
output/subway_n02_tokyo23_rough_with_shafts_round.stl
output/subway_n02_tokyo23_rough_all_lines_round.stl
output/subway_n02_tokyo23_station_shafts_round.stl
output/subway_n02_tokyo23_by_line_round/*.stl
normalized_subway_n02_tokyo23/
data/n02/N02-24_GML.zip
```

## Current result

```text
n02 lines: 603
n02 stations: 284
tube cylinders: 2963
shafts: 284
```

Route station counts:

```text
Toei Subway Asakusa Line: 20
Toei Subway Mita Line: 27
Toei Subway Oedo Line: 38
Toei Subway Shinjuku Line: 21
Tokyo Metro Chiyoda Line: 20
Tokyo Metro Fukutoshin Line: 11
Tokyo Metro Ginza Line: 19
Tokyo Metro Hanzomon Line: 14
Tokyo Metro Hibiya Line: 22
Tokyo Metro Marunouchi Line: 28
Tokyo Metro Namboku Line: 19
Tokyo Metro Tozai Line: 21
Tokyo Metro Yurakucho Line: 24
```

## Notes

- Tube and shaft radii are specified in meters and converted by `radius_m * scale`.
- `tube-sides 16` and `shaft-sides 16` make the tubes rounder than the earlier 8-sided preview.
- If `--origin-json` is not supplied, the script computes an origin from subway data bounds. This is good for subway-only preview, but not exact overlay with a building STL.
- For final overlay, pass the same origin JSON used by the building model.

## Map-aligned origin

For map-data alignment, use:

```text
origin_epsg6677_zero.json
```

This fixes the STL origin to EPSG:6677 zero:

```text
stl_x = x_epsg6677 * scale
stl_y = y_epsg6677 * scale
stl_z = z_m * scale * vertical_exaggeration
```

The Tokyo 23 wards map-origin model was generated with:

```bash
/root/.venvs/tokyo-stl-models/bin/python subway_n02_rough_model_standalone.py \
  --n02-zip data/n02/N02-24_GML.zip \
  --normalized-out-dir normalized_subway_n02_tokyo23_map_origin \
  --origin-json origin_epsg6677_zero.json \
  --out output/subway_n02_tokyo23_map_origin_with_shafts_round.stl \
  --tube-out output/subway_n02_tokyo23_map_origin_all_lines_round.stl \
  --shaft-out output/subway_n02_tokyo23_map_origin_station_shafts_round.stl \
  --by-line-dir output/subway_n02_tokyo23_map_origin_by_line_round \
  --west 139.56 \
  --east 139.93 \
  --south 35.52 \
  --north 35.82 \
  --clip-mode bbox \
  --scale 0.0005 \
  --tube-radius-m 10.0 \
  --shaft-radius-m 16.0 \
  --tube-sides 16 \
  --shaft-sides 16 \
  --min-point-spacing-m 50.0 \
  --vertical-exaggeration 2.0 \
  --station-shafts \
  --split-by-line
```

# Standalone Subway Rough Model

This note describes `subway_rough_model_standalone.py`.

The file is intentionally independent from the main CLI to reduce merge conflicts with other Codex edits.

## Why this exists

The first subway STL used `tube-radius` and `shaft-radius` directly as STL units. With `scale=0.0005`, a radius of `1.0` STL unit is visually huge compared with the model extent, so the output can look like one solid blob.

The standalone script treats radii as meters:

```text
model_radius = radius_m * scale
```

For example:

```text
tube_radius_m = 10.0
scale = 0.0005
model radius = 0.005
```

## Command used

```bash
/root/.venvs/tokyo-stl-models/bin/python subway_rough_model_standalone.py \
  --normalized-dir normalized_subway \
  --out output/subway_rough_with_shafts_mscale.stl \
  --tube-out output/subway_rough_all_lines_mscale.stl \
  --shaft-out output/subway_station_shafts_mscale.stl \
  --by-line-dir output/subway_rough_by_line_mscale \
  --scale 0.0005 \
  --tube-radius-m 10.0 \
  --shaft-radius-m 16.0 \
  --tube-sides 8 \
  --shaft-sides 8 \
  --vertical-exaggeration 2.0 \
  --station-shafts \
  --split-by-line
```

## Generated local files

These files are ignored by Git:

```text
output/subway_rough_with_shafts_mscale.stl
output/subway_rough_all_lines_mscale.stl
output/subway_station_shafts_mscale.stl
output/subway_rough_by_line_mscale/*.stl
```

## Current limitation

If `--origin-json` is not supplied, the script computes an origin from subway data bounds. This is fine for subway-only inspection, but not for exact overlay with PLATEAU building STL.

For final overlay, pass the same origin JSON used by the building model:

```bash
--origin-json origin.json
```


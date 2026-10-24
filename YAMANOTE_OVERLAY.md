# Yamanote Line Overlay

This note describes `yamanote_n02_overlay.py`.

The tube builder is shared with `subway_n02_welded_tube_model.py`; keep both
files in the same checkout. The shared implementation fixes cap winding,
parallel-transports cross sections, and closes loops with shared indices.
This does not make the estimated elevation profile surveyed data. Legacy
overlay commands still require explicit matching vertical-exaggeration and
Z-offset arguments in addition to the shared origin JSON.

The script builds a JR East Yamanote operational loop from N02 railway station data, using the same EPSG:6677 origin and scale as the subway model.

N02 does not store the full operational Yamanote loop only as `Yamanote Line`. The western arc is listed as Yamanote, while the east and southeast arcs appear under Tohoku and Tokaido railway names in the source data. The script therefore uses the 30 operational Yamanote station sequence and picks station coordinates from those three N02 source lines.

## Command Used

```bash
/root/.venvs/tokyo-stl-models/bin/python yamanote_n02_overlay.py \
  --n02-zip data/n02/N02-24_GML.zip \
  --origin-json origin_epsg6677_zero.json \
  --subway-stl output/subway_n02_tokyo23_welded_z20_r3_with_shafts_boxes_thin.stl \
  --out output/yamanote_n02_tokyo23_overlay.stl \
  --station-out output/yamanote_n02_tokyo23_stations.stl \
  --with-stations-out output/yamanote_n02_tokyo23_overlay_with_stations.stl \
  --combined-out output/subway_yamanote_n02_tokyo23_combined.stl \
  --west 139.56 \
  --east 139.93 \
  --south 35.52 \
  --north 35.82 \
  --clip-mode bbox \
  --scale 0.0005 \
  --use-estimated-elevation \
  --vertical-exaggeration 20.0 \
  --tube-radius-m 24.0 \
  --station-radius-m 45.0 \
  --tube-sides 16 \
  --min-point-spacing-m 50.0 \
  --stitch-tolerance-m 8.0 \
  --station-markers
```

## Outputs

```text
output/yamanote_n02_tokyo23_overlay.stl
output/yamanote_n02_tokyo23_stations.stl
output/yamanote_n02_tokyo23_overlay_with_stations.stl
output/subway_yamanote_n02_tokyo23_combined.stl
```

## Current Result

```text
yamanote line parts: 1
stitched paths: 1
stations: 30
z range m: 2.00 to 37.00
```

The Yamanote overlay uses an estimated station elevation profile by default. N02 has no rail height or terrain Z values, so the built-in profile is an approximate visual alignment profile, not surveyed track elevation.

With `vertical_exaggeration=20.0` and `scale=0.0005`, a 1 m elevation difference becomes `0.01` model units. The current estimated station profile ranges from 2 m to 37 m, so the Yamanote rail is visibly sloped above the underground subway model.

To provide checked elevations, pass a CSV with:

```text
station,z_m,source,note
```

and run with:

```bash
--elevation-csv yamanote_elevation.csv
```

Use `--no-use-estimated-elevation --z-m 18.0` to restore the old fixed-height behavior.

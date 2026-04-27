# Yamanote Line Overlay

This note describes `yamanote_n02_overlay.py`.

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
  --z-m 18.0 \
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
```

The Yamanote overlay uses `z_m=18.0`, so with `vertical_exaggeration=20.0` and `scale=0.0005`, the rail is placed around model Z `0.18`. The subway model remains underground with negative Z values.

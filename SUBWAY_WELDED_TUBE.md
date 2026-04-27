# Welded Subway Tube Model

This note describes `subway_n02_welded_tube_model.py`.

The earlier N02 rough model used one cylinder per line segment. That can show visible gaps at bends and joins because adjacent cylinders are not welded.

This standalone script builds each subway path as one connected tube mesh:

- Build one circular ring at each sampled line point.
- Connect neighboring rings with quad faces split into triangles.
- Cap both ends of each path.
- Stitch nearby N02 feature endpoints before tube generation.

## Command Used

```bash
/root/.venvs/tokyo-stl-models/bin/python subway_n02_welded_tube_model.py \
  --n02-zip data/n02/N02-24_GML.zip \
  --normalized-out-dir normalized_subway_n02_tokyo23_welded_z10_r3 \
  --origin-json origin_epsg6677_zero.json \
  --out output/subway_n02_tokyo23_welded_z10_r3_with_shafts.stl \
  --tube-out output/subway_n02_tokyo23_welded_z10_r3_all_lines.stl \
  --shaft-out output/subway_n02_tokyo23_welded_z10_r3_station_shafts.stl \
  --junction-box-out output/subway_n02_tokyo23_welded_z10_r3_junction_boxes.stl \
  --by-line-dir output/subway_n02_tokyo23_welded_z10_r3_by_line \
  --west 139.56 \
  --east 139.93 \
  --south 35.52 \
  --north 35.82 \
  --clip-mode bbox \
  --scale 0.0005 \
  --tube-radius-m 30.0 \
  --shaft-radius-m 48.0 \
  --junction-box-size-m 80.0 \
  --tube-sides 16 \
  --shaft-sides 16 \
  --min-point-spacing-m 50.0 \
  --stitch-tolerance-m 8.0 \
  --vertical-exaggeration 10.0 \
  --station-shafts \
  --junction-boxes \
  --split-by-line
```

## Result

```text
n02 lines: 603
stitched tube paths: 18
n02 stations: 284
tube meshes: 18
shafts: 284
junction boxes: 284
tube_radius_model: 0.015
```

Each stitched tube path was checked as watertight individually. The combined STL may still report non-watertight because multiple independently closed route meshes and station shafts intersect or overlap in one file.

Use this output when inspecting tube continuity:

```text
output/subway_n02_tokyo23_welded_z10_r3_with_shafts_boxes.stl
output/subway_n02_tokyo23_welded_z10_r3_all_lines.stl
output/subway_n02_tokyo23_welded_z10_r3_junction_boxes.stl
output/subway_n02_tokyo23_welded_z10_r3_by_line/*.stl
```

## Junction boxes

The shaft-to-tube visual intersection can still look rough because the shaft cylinders and subway tube meshes are separate solids. Use `--junction-boxes` to add one box at each station underground point.

The box is intentionally larger than the tube/shaft contact area. It is a visual connector, not an actual underground structure.

For shaft radius 48 m, prefer oriented boxes larger than the shaft diameter:

```bash
--junction-box-size-m 0 \
--junction-box-length-m 220.0 \
--junction-box-width-m 150.0 \
--junction-box-height-m 120.0
```

When `--junction-box-size-m 0` is used, the script creates a rectangular box and rotates it in the horizontal plane to match the nearest subway tube segment.

Current oriented-box output:

```text
output/subway_n02_tokyo23_welded_z10_r3_with_shafts_boxes_oriented.stl
output/subway_n02_tokyo23_welded_z10_r3_junction_boxes_oriented.stl
```

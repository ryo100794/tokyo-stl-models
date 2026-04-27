# Geology Layer STL

This note describes `geology_layer_model.py`.

The script creates separate geology STL files that align with the subway model by using the same EPSG:6677 origin, scale, bbox, and vertical exaggeration.

The current default data is an approximate Tokyo 23 wards profile. It is intended only as a visual placeholder until checked boundary surfaces from GSJ Urban Geological Map or another source are converted into CSV.

## Command Used

```bash
/root/.venvs/tokyo-stl-models/bin/python geology_layer_model.py \
  --origin-json origin_epsg6677_zero.json \
  --out-dir output/geology_tokyo23_z50 \
  --west 139.56 \
  --east 139.93 \
  --south 35.52 \
  --north 35.82 \
  --grid-x 64 \
  --grid-y 52 \
  --scale 0.0005 \
  --vertical-exaggeration 50.0 \
  --explode-offset-model 0.35
```

## Outputs

```text
output/geology_tokyo23_z50/geology_boundaries_combined.stl
output/geology_tokyo23_z50/geology_layers_combined.stl
output/geology_tokyo23_z50/geology_boundaries_exploded.stl
output/geology_tokyo23_z50/geology_layers_exploded.stl
output/geology_tokyo23_z50/geology_boundary_markers_exploded.stl
output/geology_tokyo23_z50/geology_boundary_surface.stl
output/geology_tokyo23_z50/geology_boundary_base_fill.stl
output/geology_tokyo23_z50/geology_boundary_base_alluvium.stl
output/geology_tokyo23_z50/geology_boundary_base_tokyo_upper.stl
output/geology_tokyo23_z50/geology_boundary_base_tokyo_lower.stl
output/geology_tokyo23_z50/geology_boundary_base_kazusa_top.stl
output/geology_tokyo23_z50/geology_layer_fill.stl
output/geology_tokyo23_z50/geology_layer_alluvium.stl
output/geology_tokyo23_z50/geology_layer_tokyo_upper.stl
output/geology_tokyo23_z50/geology_layer_tokyo_lower.stl
output/geology_tokyo23_z50/geology_layer_kazusa_top.stl
output/geology_tokyo23_z50/geology_boundary_samples.csv
```

## CSV Boundary Input

To replace the approximate profile with checked boundary points, pass:

```bash
--surface-csv geology_boundary_points.csv
```

CSV columns:

```text
boundary,lon,lat,z_m
```

Supported boundary names:

```text
surface
base_fill
base_alluvium
base_tokyo_upper
base_tokyo_lower
base_kazusa_top
```

The script uses inverse-distance interpolation from CSV points and falls back to the built-in rough profile where points are missing.

## Notes

The default model is not a surveyed geological model. It is a separate STL scaffold for checking how geological boundary surfaces and layer volumes can be displayed with the subway model.

Use GSJ Urban Geological Map boundary surfaces or checked borehole-derived data before treating the shape as geology.

If the combined layer STL looks like one solid block in a single-color viewer, use:

```text
geology_layers_exploded.stl
geology_boundaries_exploded.stl
geology_boundary_markers_exploded.stl
```

These files offset each boundary or layer diagonally in model space so the layer order can be inspected. They are for checking the layer stack, not for direct geographic overlay.

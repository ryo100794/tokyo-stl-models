# Codex Handoff

This repository contains a Python CLI for converting Project PLATEAU CityGML building data to STL.

## Repository

- GitHub: https://github.com/ryo100794/tokyo-stl-models
- Default branch: main
- Main script: `plateau_citygml_to_stl.py`
- README: `README.md`

## Important data policy

Do not commit downloaded PLATEAU data or generated STL files.

The repository intentionally tracks only:

- Source code
- README and handoff docs
- `data/.gitkeep`
- `output/.gitkeep`

The following are ignored:

- `*.zip`
- `*.stl`
- `data/*`
- `output/*`

## Local data already downloaded

On the working machine used for this handoff, bbox-targeted building GML files were downloaded under:

```text
data/tokyo23_actual/udx/bldg/
```

There are 88 `.gml` files, about 5.0 GB total.

The URL manifest is:

```text
data/actual_gml_urls.json
```

This file is ignored and is not pushed.

The target bbox used so far is:

```text
west  139.655
east  139.735
south 35.640
north 35.715
```

This covers the Shinjuku, Shibuya, and Shimokitazawa area used in the original examples.

## Local STL outputs already generated

These outputs exist locally but are ignored by Git:

```text
output/test_lod1_5files.stl
output/lod1_batch_0005.stl
output/lod1_batch_0010.stl
output/lod1_batch_0015.stl
```

They were generated with `--lod 1` and small file batches.

Important: these were generated with `--base-thickness 3.0`, so each file contains its own flat base plate. They are useful for checking geometry, but they are not ideal final merge inputs.

For final production, regenerate parts with:

```text
--base-thickness 0
```

Then merge building-only STL parts and add one final base plate.

## Memory situation

The machine had about:

```text
Mem total: 10 GiB
Swap total: 6 GiB
```

During CityGML processing, swap became nearly full. Do not run all 88 files at once.

Safe operating pattern:

1. Check memory before each batch:

```bash
free -h
```

2. Use LOD1 and small file batches:

```bash
/root/.venvs/tokyo-stl-models/bin/python plateau_citygml_to_stl.py stl \
  --root data/tokyo23_actual \
  --out output/parts/lod1_0018.stl \
  --west 139.655 \
  --east 139.735 \
  --south 35.640 \
  --north 35.715 \
  --scale 0.0005 \
  --base-thickness 0 \
  --lod 1 \
  --file-offset 18 \
  --max-files 2
```

3. Check memory again after each batch:

```bash
free -h
```

If swap free is close to zero, pause.

## Recommended next steps

1. Add `merge-stl` command.

Use streaming binary STL merging:

- Read 80-byte header and 4-byte face count from each binary STL.
- Sum face counts.
- Write one output header and total face count.
- Copy triangle records in chunks.

This avoids loading all meshes into memory.

2. Add `add-base` command.

For a final merged STL:

- Stream-read all triangles to compute bounds.
- Stream-copy triangles to a new STL.
- Append one rectangular base plate around the final bounds.

3. Regenerate parts without base plates.

Use `--base-thickness 0`.

Suggested chunking on this machine:

```text
--max-files 2
```

Start from offset 0 if making clean final parts. The previous local outputs with base plates should not be used for the final merged model.

4. Merge clean parts.

After `merge-stl` is implemented, merge the generated `output/parts/*.stl`.

5. Add a single base plate.

After `add-base` is implemented, add one base plate to the final merged STL.

## Current CLI notes

The current script supports:

- `download`
- `unzip`
- `stl`

The `stl` command supports low-memory operation by default. It scans once to compute bounds and face counts, then writes binary STL in a second pass.

Useful options:

```text
--lod 1
--file-offset N
--max-files N
--skip-files-larger-than-mb N
--base-thickness 0
```

Avoid `--repair` for large runs because it forces in-memory mesh processing through trimesh.

## Validation already run

The script passed:

```bash
/root/.venvs/tokyo-stl-models/bin/python -m py_compile plateau_citygml_to_stl.py
```

Small synthetic CityGML STL output also passed earlier.

## Notes about boundaries

PLATEAU CityGML files are usually split by standard regional mesh tiles, not by exact ward boundary lines.

When combining ward data:

- Deduplicate by mesh code or URL.
- Use the same bbox, scale, and origin strategy for every batch.
- Avoid per-batch base plates.

## Notes about ground height

Building CityGML contains Z values for building geometry, including surfaces near ground level, but it is not a continuous terrain mesh.

The current base plate is flat. For terrain relief, add support for terrain/DEM data separately.


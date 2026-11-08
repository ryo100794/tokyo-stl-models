#!/usr/bin/env python3
"""Reproduce R14 rails in a fresh directory from saved inputs and bounded repairs.

Generated geometry is a candidate until its own serialized bytes pass the
appropriate independent checks. This command never declares a print release.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
LINE_IDS = [f'L{i:02d}' for i in range(1, 15)]
PATCHED_IDS = {'L07', 'L09', 'L12', 'L13', 'L14'}
SOURCE_FILES = [
    'reproduce_rails_r14.py', 'build_r14.py', 'geometry_r14.py',
    'rail_components_r14.py', 'export_repaired_parts.py', 'robust_export_r14.py',
    'round_internal_junctions_r14.py', 'validate_saved_meshes.py',
    'verify_saved_wall_regions.py', 'verify_installed_channels.py',
    'mesh_ray_checks_r13.py', 'check_stored_stl_normals.py',
    'inspect_print_orientation.py', 'model_coordinates.py',
    'render_r14.py', 'cpu_mesh_render.py', 'cpu_mesh_raster.cpp',
    'build_flow_coupon_R14.py', 'requirements_rebuild_R14.txt',
]
BODY_SOURCE_FILES = ['reproduce_body_r14.py', 'extract_principal_boundary.py', 'repair_selected_seams.py', 'core_signature_r14.py']
PACKAGES = ['numpy', 'scipy', 'trimesh', 'manifold3d', 'pymeshlab', 'vtk',
            'rtree', 'networkx', 'Pillow', 'lxml', 'requests', 'pyproj', 'PyYAML']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def dump(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')


def file_record(path, root):
    return {'path': path.relative_to(root).as_posix(),
            'bytes': path.stat().st_size, 'sha256': sha(path)}


def require_parts_mesh_checks(report, models):
    expected = {p.name: sha(p) for p in models}
    tests = report.get('tests', [])
    if len(tests) != len(expected) or {t.get('basename') for t in tests} != set(expected):
        raise ValueError('Independent validation does not cover exactly the current 15 STL files')
    failures = []
    for row in tests:
        name = row['basename']
        good = row.get('sha256') == expected[name]
        good &= row.get('finite_vertex_coordinates') is True
        good &= row.get('trimesh_is_watertight') is True
        good &= row.get('trimesh_is_winding_consistent') is True
        good &= row.get('signed_volume_mm3', 0) > 0
        for field in ['boundary_edges', 'edges_with_more_than_two_faces',
                      'inconsistent_two_face_edge_orientation', 'zero_area_triangles_exact',
                      'duplicate_faces_ignoring_orientation']:
            good &= row.get(field) == 0
        good &= row.get('vertex_links', {}).get('nonmanifold_vertex_count') == 0
        good &= row.get('self_intersection_test', {}).get('selected_face_count') == 0
        if not good:
            failures.append(name)
    if failures:
        raise ValueError(f'Saved part mesh checks failed or were incomplete: {failures}')


def integrate_cusps(root, hull_cache, cusp_out, final_cache):
    report = read(cusp_out / 'internal_cusp_repair_trials.json')
    manifest = read(cusp_out / 'artifact_manifest_R14_internal_junctions.json')
    rows = report.get('lines', [])
    if {row['line_id'] for row in rows} != PATCHED_IDS:
        raise ValueError('The saved R14 repair contract requires all five specified line repairs')
    if report.get('numerical_integration_gates_pass') is not True:
        raise ValueError('Cusp integration numerical gates failed or are missing')
    if manifest.get('numerical_integration_gates_pass') is not True:
        raise ValueError('Cusp artifact manifest has no passing integration gate')
    if len(report.get('all_91_final_core_pairs', [])) != 91:
        raise ValueError('All 91 final inter-line core pairs must be checked')
    if report.get('all_final_core_pair_gaps_ge_0_5_mm') is not True:
        raise ValueError('A final inter-line core separation guard failed')
    listed = {}
    for entry in manifest['files']:
        path = (cusp_out / entry['path']).resolve()
        if not path.is_relative_to(cusp_out.resolve()) or not path.is_file():
            raise ValueError(f'Invalid cusp manifest path: {entry["path"]}')
        if path.stat().st_size != entry['bytes'] or sha(path) != entry['sha256']:
            raise ValueError(f'Cusp artifact hash mismatch: {entry["path"]}')
        listed[entry['path']] = entry
    for row in rows:
        lid = row['line_id']
        if row.get('numerical_integration_gates_pass') is not True:
            raise ValueError(f'Line integration guard failed: {lid}')
        for kind in ['core', 'outer']:
            if sha(hull_cache / f'{lid}_{kind}.npz') != row[f'original_{kind}_sha256']:
                raise ValueError(f'Line input changed since the cusp calculation: {lid}/{kind}')
        candidate = cusp_out / f'rail_R14_{lid}_cusp_trial.stl'
        if candidate.name not in listed or sha(candidate) != row['saved_stl']['sha256']:
            raise ValueError(f'Line candidate not covered by validated hashes: {lid}')
        core = f'{lid}_core_after.npz'
        if core not in listed:
            raise ValueError(f'Final core not covered by the artifact manifest: {lid}')
    final_cache.mkdir(parents=True, exist_ok=False)
    for lid in LINE_IDS:
        source = cusp_out / f'{lid}_core_after.npz' if lid in PATCHED_IDS else hull_cache / f'{lid}_core.npz'
        shutil.copy2(source, final_cache / f'{lid}_core.npz')
    # Geometry imports resolve to the fixed source snapshot copied into ROOT/source.
    sys.path.insert(0, str(root / 'source'))
    from geometry_r14 import load_npz, save_npz, unite
    all_cores = unite([load_npz(final_cache / f'{lid}_core.npz') for lid in LINE_IDS])
    save_npz(all_cores, final_cache / 'all_cores.npz')
    for lid in sorted(PATCHED_IDS):
        shutil.copy2(cusp_out / f'rail_R14_{lid}_cusp_trial.stl', root / 'models' / f'rail_R14_{lid}.stl')
    receipt = {
        'phase': 'numerically guarded individual-line integration; no print release',
        'source_report_sha256': sha(cusp_out / 'internal_cusp_repair_trials.json'),
        'source_manifest_sha256': sha(cusp_out / 'artifact_manifest_R14_internal_junctions.json'),
        'replaced_lines': sorted(PATCHED_IDS),
        'unchanged_core_lines': [lid for lid in LINE_IDS if lid not in PATCHED_IDS],
        'cores': [file_record(final_cache / f'{lid}_core.npz', root) for lid in LINE_IDS],
        'all_cores': file_record(final_cache / 'all_cores.npz', root),
        'models': [file_record(root / 'models' / f'rail_R14_{lid}.stl', root) for lid in LINE_IDS],
        'global_wall_thickness_proven': False,
        'whole_body_acceptance': 'pending',
        'physical_tests_performed': False,
    }
    dump(root / 'reports/final_line_integration_R14.json', receipt)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', type=Path, required=True,
                   help='Read-only saved input root containing data/, reference/ and the unchanged coupon')
    p.add_argument('--out', type=Path, required=True, help='New empty working directory')
    p.add_argument('--include-body', action='store_true', help='Also generate and inspect the whole-body candidate')
    p.add_argument('--rebuild-positive', action='store_true',
                   help='Replay positive-material seed repairs instead of the pinned positive checkpoint')
    p.add_argument('--render', action='store_true', help='Create optional CPU previews; requires C++17 and DejaVu Sans')
    a = p.parse_args()
    if a.rebuild_positive and not a.include_body:
        p.error('--rebuild-positive requires --include-body')
    source = Path(__file__).resolve().parent
    inp, root = a.input_root.resolve(), a.out.resolve()
    if inp == root or root.is_relative_to(inp) or inp.is_relative_to(root):
        p.error('Input and output roots must be separate, non-nested directories')
    if root.exists() and any(root.iterdir()):
        p.error('Output must be new or empty; existing results are never overwritten')
    required = ['data/rail_paths_R14.json', 'data/ports_R14.json',
                'data/coordinate_systems_R14.json', 'data/internal_cusp_seed_R14.json',
                'data/flow_coupon_R14_specification.json', 'models/flow_coupon_R14.stl',
                'reference/subway_R12_lite.glb']
    if a.include_body:
        required.append('data/body_repair_recipe_R14.json')
        required.append('data/body_core_geometry_R14.json')
    for relative in required:
        if not (inp / relative).is_file():
            p.error(f'Missing saved input: {relative}')
    runtime = SOURCE_FILES + (BODY_SOURCE_FILES if a.include_body else [])
    for name in runtime:
        if not (source / name).is_file():
            p.error(f'Missing current source dependency: {name}')
    root.mkdir(parents=True, exist_ok=True)
    shutil.copytree(inp / 'data', root / 'data')
    shutil.copytree(inp / 'reference', root / 'reference')
    (root / 'source').mkdir(); (root / 'models').mkdir(); (root / 'reports').mkdir()
    for name in runtime:
        shutil.copy2(source / name, root / 'source' / name)
    shutil.copy2(inp / 'models/flow_coupon_R14.stl', root / 'models/flow_coupon_R14.stl')
    started = time.time()
    status = {
        'schema': 'R14-reproduction-status-v1',
        'phase': 'R14 regeneration candidates; independent whole-body acceptance required; NOT a print release',
        'input_root_read_only': str(inp), 'output_root': str(root),
        'inputs': [file_record(path, root) for folder in ['data', 'reference']
                   for path in sorted((root / folder).rglob('*')) if path.is_file()],
        'sources': [file_record(root / 'source' / name, root) for name in runtime],
        'coupon_copied_unchanged': file_record(root / 'models/flow_coupon_R14.stl', root),
        'runtime': {'python': sys.version, 'packages': {name: importlib.metadata.version(name) for name in PACKAGES}},
        'commands': [], 'individual_mesh_checks': 'pending', 'whole_body_acceptance': 'pending',
        'global_wall_thickness_proven': False, 'physical_tests_performed': False,
    }
    def persist():
        status['seconds'] = time.time() - started
        dump(root / 'reports/pipeline_status_R14.json', status)
    def run(script, *arguments):
        command = [sys.executable, str(root / 'source' / script), *map(str, arguments)]
        row = {'argv': command, 'status': 'running'}; status['commands'].append(row); persist()
        code = subprocess.call(command, cwd=root)
        row.update(returncode=code, status='completed' if code == 0 else 'failed'); persist()
        if code:
            raise RuntimeError(f'{script} exited {code}; inspect its report and retained candidates')
    hull, final, cusp = root / 'cache/hull', root / 'cache/final', root / 'reports/internal_junctions'
    try:
        run('build_r14.py', '--root', root, '--cache', hull, '--stage', 'rails')
        run('build_r14.py', '--root', root, '--cache', hull, '--stage', 'parts')
        # This report describes only the baseline parts before five guarded replacements.
        (root / 'reports/repaired_part_exports_R14.json').rename(root / 'reports/baseline_part_exports_R14.json')
        run('round_internal_junctions_r14.py', '--root', root, '--cache', hull, '--out', cusp)
        integrate_cusps(root, hull, cusp, final)
        parts = [root / 'models/flow_coupon_R14.stl'] + [root / 'models' / f'rail_R14_{lid}.stl' for lid in LINE_IDS]
        run('validate_saved_meshes.py', *parts, '--out', root / 'reports/independent_parts')
        require_parts_mesh_checks(read(root / 'reports/independent_parts/independent_saved_meshes_report.json'), parts)
        status['individual_mesh_checks'] = 'passed_for_current_saved_bytes_with_recorded_predicate_limits'
        run('verify_saved_wall_regions.py', '--root', root, '--out', root / 'reports/saved_wall_regions_R14.json')
        wall = read(root / 'reports/saved_wall_regions_R14.json')
        if len(wall.get('lines', [])) != 14 or {row.get('line_id') for row in wall['lines']} != set(LINE_IDS):
            raise ValueError('Finite wall screening did not cover all 14 current line files')
        wall_clear = all(row.get('sha256') == sha(root / 'models' / f'rail_R14_{row["line_id"]}.stl')
                         and all(row.get(k) == 0 for k in ['radial_below_threshold_count', 'normal_below_threshold_count', 'missing_boundary'])
                         for row in wall['lines'])
        status['finite_ordinary_outer_wall_screen'] = 'clear_with_sampling_limits' if wall_clear else 'requires_review'
        run('check_stored_stl_normals.py', *parts, '--out', root / 'reports/stored_stl_normals_R14.json')
        run('inspect_print_orientation.py', '--root', root, '--out', root / 'reports/orientation')
        if a.include_body:
            if a.rebuild_positive:
                run('build_r14.py', '--root', root, '--cache', hull, '--stage', 'base')
                run('build_r14.py', '--root', root, '--cache', hull, '--stage', 'assemble')
            body_out = root / 'body_reproduction'
            arguments = ['--root', root, '--cache', hull, '--core-cache', final, '--out', body_out]
            if a.rebuild_positive:
                arguments.append('--rebuild-positive')
            run('reproduce_body_r14.py', *arguments)
            record = read(body_out / 'body_reproduction_R14.json')
            body = body_out / 'subway_R14_hollow_candidate_print.stl'
            if record.get('candidate_sha256') != sha(body):
                raise ValueError('Body candidate bytes do not match its reproduction report')
            status['body_candidate'] = file_record(body, root)
            run('validate_saved_meshes.py', body, '--out', root / 'reports/independent_body')
            run('verify_installed_channels.py', '--root', root, '--model', body,
                '--out', root / 'reports/installed_channels_R14.json')
            run('check_stored_stl_normals.py', body, '--out', root / 'reports/body_stored_normals_R14.json')
            status['whole_body_acceptance'] = 'pending exact all-triangle predicate and removed-surface containment review for this SHA256'
            status['whole_body_report_bundle'] = [file_record(path, root) for folder in ['body_reproduction', 'reports/independent_body']
                                                   for path in sorted((root / folder).glob('*.json'))]
        if a.render:
            run('render_r14.py', '--root', root)
        status['pipeline_execution'] = 'completed; numerical results retain their stated scope'
    except Exception as error:
        status['pipeline_execution'] = 'failed_or_incomplete'; status['error'] = repr(error)
        persist(); raise
    persist()
    print(json.dumps({'status_report': str(root / 'reports/pipeline_status_R14.json'),
                      'whole_body_acceptance': status['whole_body_acceptance'],
                      'physical_tests_performed': False}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()

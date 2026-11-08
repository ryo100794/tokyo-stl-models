#!/usr/bin/env python3
"""Build a one-piece flow-and-colour-isolation coupon for rail holder R14.

Run with the R14 source directory and installed dependencies on PYTHONPATH.
The coupon is a digital design. It has not been printed or physically tested.
All positive solids, including the substrate, support pedestals and labels,
are Boolean-unioned before ALL lumen cutters are subtracted. This avoids a
closed substrate or label silently obstructing a nominally hollow tube.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import manifold3d as mf
import trimesh

from geometry_r14 import capsule_path, cylinder_between, unite, arrays, write_stl


SIDES = 24
PLATE = [54., 68., 1.2]
ROUTE_Z = 2.2
WALL = .6
PORT_LENGTH = 3.6
PORT_RI = 1.
PORT_RO = 1.6


def plate():
    """Make a shallow rounded rectangular substrate with 2 mm corner radii."""
    width, height, thickness = PLATE
    radius = 2.
    pieces = [mf.Manifold.cube([width - 2 * radius, height, thickness]).translate([radius, 0, 0]),
              mf.Manifold.cube([width, height - 2 * radius, thickness]).translate([0, radius, 0])]
    pieces.extend(mf.Manifold.cylinder(thickness, radius, radius, SIDES).translate([x, y, 0])
                  for x in (radius, width - radius) for y in (radius, height - radius))
    return unite(pieces)


def u_path(y0):
    """An exactly 80 mm saved polyline with a 3 mm-radius semicircular bend."""
    angles = np.linspace(-math.pi / 2, math.pi / 2, 49)
    arc = np.column_stack([3. * np.cos(angles), 3. * np.sin(angles), np.zeros(len(angles))])
    arc_length = float(np.linalg.norm(np.diff(arc, axis=0), axis=1).sum())
    run = (80. - arc_length) / 2
    cx, cy = 8. + run, y0 + 3.
    arc += [cx, cy, ROUTE_Z]
    return np.vstack([[8., y0, ROUTE_Z], arc, [8., y0 + 6., ROUTE_Z]])


def nozzle(point, radius):
    """Match the main model's 2 mm neck and linear taper to its common mouth."""
    point = np.asarray(point, float)
    axis = np.array([0., 0., 1.])
    a, mouth = point + 2. * axis, point + PORT_LENGTH * axis
    outer = unite([cylinder_between(point - .15 * axis, a, radius + WALL, sides=SIDES),
                   cylinder_between(a, mouth, radius + WALL, PORT_RO, sides=SIDES)])
    inner = unite([cylinder_between(point - .20 * axis, a, radius, sides=SIDES),
                   cylinder_between(a, mouth, radius, PORT_RI, sides=SIDES),
                   cylinder_between(mouth, mouth + .8 * axis, PORT_RI, sides=SIDES),
                   mf.Manifold.sphere(radius, SIDES).translate(point)])
    return outer, inner, mouth


DIGITS = {
    '0': 'ab cdef'.replace(' ', ''), '1': 'bc', '2': 'abged',
    '3': 'abgcd', '4': 'fgbc', '5': 'afgcd', '6': 'afgecd',
    '7': 'abc', '8': 'abcdefg', '9': 'abfgcd', 'A': 'abcefg', 'b': 'cdefg',
}


def raised_label(text, origin):
    """Use reproducible, legible 3.2 mm-high geometric digits and letters."""
    x, y = origin
    w, h, stroke, depth, baseline = 1.8, 3.2, .32, .34, PLATE[2] - .04
    pieces = []
    for char in text:
        if char == '.':
            pieces.append(mf.Manifold.cylinder(depth, .18, .18, 20).translate([x + .20, y + .18, baseline]))
            x += .8
            continue
        segments = {
            'a': (x, y + h - stroke, w, stroke),
            'g': (x, y + (h - stroke) / 2, w, stroke),
            'd': (x, y, w, stroke),
            'f': (x, y + h / 2, stroke, h / 2),
            'b': (x + w - stroke, y + h / 2, stroke, h / 2),
            'e': (x, y, stroke, h / 2),
            'c': (x + w - stroke, y, stroke, h / 2),
        }
        for name in DIGITS[char]:
            sx, sy, sw, sh = segments[name]
            pieces.append(mf.Manifold.cube([sw, sh, depth]).translate([sx, sy, baseline]))
        x += w + .65
    return unite(pieces)


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    positives = [plate()]
    core_solids, records = [], []
    channel_specs = [(f'D{int(d * 10):02d}', d, u_path(y), f'{d:.1f}')
                     for d, y in zip([1.0, 1.2, 1.4], [6., 18., 30.])]
    channel_specs += [
        ('cross_A', 1.2, np.array([[8., 53., ROUTE_Z], [46., 53., ROUTE_Z]]), 'A'),
        ('cross_B', 1.2, np.array([[27., 44., ROUTE_Z + 1.9], [27., 62., ROUTE_Z + 1.9]]), 'b'),
    ]
    for index, (identity, diameter, points, label) in enumerate(channel_specs):
        print('Building channel', identity, flush=True)
        radius = diameter / 2
        positives.append(capsule_path(points, radius + WALL, sides=SIDES, tolerance=0.))
        channel_core = [capsule_path(points, radius, sides=SIDES, tolerance=0.)]
        ports = []
        for end, point in [('start', points[0]), ('end', points[-1])]:
            outer, inner, mouth = nozzle(point, radius)
            positives.append(outer)
            channel_core.append(inner)
            ports.append({'end': end, 'attachment_mm': point.tolist(), 'mouth_center_mm': mouth.tolist(),
                          'axis': [0., 0., 1.], 'length_mm': PORT_LENGTH,
                          'inner_mouth_diameter_mm': 2., 'outer_mouth_diameter_mm': 3.2,
                          'straight_neck_length_mm': 2., 'linear_taper_length_mm': 1.6})
            if identity == 'cross_B':
                positives.append(cylinder_between([point[0], point[1], 1.], point, 1.35, sides=SIDES))
        core = unite(channel_core)
        core_solids.append(core)
        if index < 3:
            positives.append(raised_label(label, [20., points[0, 1] + 1.4]))
        else:
            positives.append(raised_label(label, [3., 49.] if identity == 'cross_A' else [30., 60.]))
        records.append({'id': identity, 'raised_marking': label,
                        'inner_diameter_mm': diameter, 'outer_diameter_mm': diameter + 1.2,
                        'nominal_radial_wall_mm': WALL,
                        'centerline_mm': points.tolist(),
                        'centerline_length_excluding_ports_mm': float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum()),
                        'ports': ports, 'core_connected_component_count': len(core.decompose())})

    # A/B two-colour crossing: the core radii are 0.6 mm and centres are 1.9 mm apart.
    # The analytical centre distance is a conservative lower bound for faceted cores.
    isolation_gap = 1.9 - 1.2
    intersection_rows = []
    for i, a in enumerate(core_solids):
        for j, b in enumerate(core_solids[i + 1:], i + 1):
            intersect = a ^ b
            row = {'channel_a': records[i]['id'], 'channel_b': records[j]['id'],
                   'intersection_volume_mm3': float(intersect.volume()),
                   'intersection_empty': bool(intersect.is_empty())}
            assert row['intersection_empty'], row
            intersection_rows.append(row)

    # Deliberately keep the order explicit: union every positive, then cut every core.
    solid_positive = unite(positives)
    solid_cores = unite(core_solids)
    print('Subtracting all five cores from the united coupon', flush=True)
    final = solid_positive - solid_cores
    assert final.status() == mf.Error.NoError
    assert len(final.decompose()) == 1
    assert len(solid_cores.decompose()) == 5
    residual = final ^ solid_cores
    residual_volume = abs(float(residual.volume()))
    print('Post-difference intersection', residual.is_empty(), residual_volume, 'mm3', flush=True)
    # Re-intersection of the two float32 surfaces can retain sub-voxel slivers
    # at the same boundary. Record their volume instead of calling them empty.
    assert residual_volume < 1e-5
    stl_path = out / 'flow_coupon_R14.stl'
    # Binary STL stores float32 coordinates. Collapse submicrometre redundant
    # Boolean seams before this rounding so they do not become zero-area faces.
    final = final.simplify(1e-5)
    qa = write_stl(final, stl_path)

    # Reload STL independently with trimesh; count the boundary of each 2 mm opening.
    saved = trimesh.load_mesh(stl_path, process=True)
    raw_saved_face_count = len(saved.faces)
    saved.update_faces(saved.nondegenerate_faces())
    saved.update_faces(saved.unique_faces())
    saved.remove_unreferenced_vertices()
    assert saved.is_watertight and saved.is_winding_consistent
    assert len(saved.split(only_watertight=False)) == 1
    stl_path.write_bytes(trimesh.exchange.stl.export_stl(saved))
    saved = trimesh.load_mesh(stl_path, process=True)
    assert saved.is_watertight and saved.is_winding_consistent
    qa.update({'reloaded_stl_watertight': bool(saved.is_watertight),
               'pre_stl_simplification_tolerance_mm': 1e-5,
               'post_stl_removed_degenerate_or_duplicate_face_count': raw_saved_face_count - len(saved.faces),
               'final_saved_triangle_count': len(saved.faces),
               'reloaded_stl_winding_consistent': bool(saved.is_winding_consistent),
               'reloaded_stl_connected_part_count': len(saved.split(only_watertight=False)),
               'all_5_core_components_separate': True,
               'core_pair_boolean_checks': intersection_rows,
               'final_intersection_with_all_cores_empty': bool(residual.is_empty()),
               'post_difference_roundoff_residual_volume_mm3': residual_volume,
               'allowed_post_difference_numeric_residual_volume_mm3': 1e-5,
               'physical_print_or_flow_test_performed': False,
               'sha256': hashlib.sha256(stl_path.read_bytes()).hexdigest()})
    spec = {
        'revision': 'R14', 'artifact': 'one_piece_flow_and_colour_isolation_coupon',
        'units': 'mm', 'coordinate_system': 'X right, Y longitudinal, Z above substrate',
        'physical_test_status': 'Not printed and not physically flow-tested; digital geometry checks only.',
        'plate_size_mm': PLATE, 'plate_corner_radius_mm': 2.,
        'overall_bounds_mm': saved.bounds.tolist(), 'overall_size_mm': saved.extents.tolist(),
        'print_part_count': 1, 'independent_channel_count': 5, 'port_count': 10,
        'channel_sections_polygon_sides': SIDES,
        'nominal_wall_definition': 'Radial CAD difference. Polygon facets reduce straight-tube wall normal distance by cos(pi/24), matching the R14 main model.',
        'minimum_nominal_straight_facet_wall_mm': WALL * math.cos(math.pi / SIDES),
        'channels': records,
        'crossing': {'channel_lower': 'cross_A', 'channel_upper': 'cross_B',
                     'crossing_xy_mm': [27., 53.], 'center_vertical_separation_mm': 1.9,
                     'analytical_minimum_lumen_to_lumen_gap_mm': isolation_gap,
                     'required_minimum_gap_mm': .6, 'meets_gap_by_construction': isolation_gap >= .6},
        'recommended_test_sequence': [
            'Print in the same material, exposure and layer settings intended for the full model.',
            'Wash all five channels through their two open ends before final post-curing, following the resin and printer cleaning instructions.',
            'Confirm wash liquid exits the opposite end of each marked 80 mm channel; record differences between 1.0, 1.2 and 1.4 mm bores.',
            'After compatible cleaning and post-curing, introduce two contrasting liquids separately through A and b; check for cross-colour contamination or leakage.',
            'Use gentle syringe pressure and material-compatible fluid; a successful digital check is not a physical leak or flow qualification.',
        ],
        'boolean_order': 'UNION(all positive substrate, tube, pedestal, nozzle and marking solids) MINUS UNION(all five full-length channel and nozzle cores)',
        'geometry_qa': qa,
    }
    (out / 'flow_coupon_R14_specification.json').write_text(json.dumps(spec, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (out / 'README_flow_coupon_R14_JA.txt').write_text(
        '地下鉄ペン立て R14 通液・別色隔離試験片\n\n'
        '台座54×68×1.2 mm、印刷部品1点、独立流路5本、開口10個です。\n'
        '刻印1.0/1.2/1.4は流路内径mmを表し、それぞれ曲がりを含む中心線長80 mmです。\n'
        '各流路の名目壁厚は0.6 mm、外径は2.2/2.4/2.6 mmです。\n'
        '全開口は内径2.0 mm、外径3.2 mm、長さ3.6 mmです。2.0 mmの直線部から1.6 mmで直線的に広がります。\n'
        'Aとbは別色確認用の上下交差流路です。内径は両方1.2 mm、交差部の内腔間隔は最低0.7 mmです。\n'
        '本体と同じ樹脂・積層厚・露光条件で印刷し、二次硬化の前に各入口から洗浄液が反対の口へ抜けるか確認してください。\n'
        '洗浄と後硬化は使用する樹脂・機器の指示に従ってください。Aとbには異なる色を別々に注入して混色と漏れを確認します。\n'
        '形状は全ての外形を一体化した後で全流路を引き抜いた実中空STLです。デジタル検査のみを実施しました。\n'
        '実際の造形、通液、耐圧、漏れの試験は未実施です。\n', encoding='utf-8')
    print(json.dumps({'stl': str(stl_path), 'qa': qa, 'size_mm': saved.extents.tolist()}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent)
    build(parser.parse_args().out)

"""Regression tests for a fixed model origin across independent CityGML tiles."""
import argparse
import contextlib
import io
import json
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from model_coordinates import ModelCoordinates
import plateau_citygml_to_stl as plateau


def stl_triangles(path):
    raw = Path(path).read_bytes()
    count, = struct.unpack_from('<I', raw, 80)
    if len(raw) != 84 + 50 * count:
        raise AssertionError('STL count differs from saved records')
    return np.array([struct.unpack_from('<12fH', raw, 84 + i * 50)[3:12]
                     for i in range(count)]).reshape(-1, 3, 3)


class ModelCoordinatesTests(unittest.TestCase):
    def test_explicit_units_and_vertical_scale(self):
        c = ModelCoordinates(-6155.9956, -34454.12365, mm_per_m=1 / 120,
                             vertical_exaggeration=200 / 3)
        np.testing.assert_allclose(c.point((-6035.9956, -34214.12365, 1.8)), [1, 2, 1], atol=1e-12)

    def test_roundtrip_and_vector_transform(self):
        c = ModelCoordinates(13, -24, 1.2, .008, 66, 4)
        source = np.array([[13, -24, 1.2], [-4500, -34000, -47.2], [0, 0, 0]])
        model = c.points(source)
        np.testing.assert_allclose(model, [c.point(point) for point in source], atol=1e-13)
        np.testing.assert_allclose([c.inverse(point) for point in model], source, atol=1e-10)

    def test_legacy_origin_and_roundtrip_json(self):
        legacy = ModelCoordinates.from_dict({'origin_x': 1, 'origin_y': 2, 'scale': .25})
        self.assertEqual(legacy.point((5, 6, 4)), (1, 1, 1))
        self.assertEqual(ModelCoordinates.from_dict(legacy.to_dict()), legacy)

    def test_contradictory_or_invalid_units_rejected(self):
        for update in ({'scale': .1, 'mm_per_m': .2}, {'model_units': 'm'},
                       {'mm_per_m': 0}, {'origin_x': float('nan')},
                       {'projected_axis_order': ['northing', 'easting', 'height']}):
            with self.subTest(update=update), self.assertRaises(ValueError):
                ModelCoordinates.from_dict({'origin_x': 1, 'origin_y': 2, **update})

    def test_declared_epsg6677_native_axes_are_normalized(self):
        transforms = plateau.make_transformers()
        ring = plateau.parse_poslist_text('-34454 -6155 3 -34454 -6143 3 -34442 -6155 6',
                                         '3', transforms, 'urn:ogc:def:crs:EPSG::6677')
        self.assertEqual(ring.points_m[0], (-6155, -34454, 3))
        undeclared = plateau.parse_poslist_text('-6155 -34454 3 -6143 -34454 3 -6155 -34442 6',
                                              '3', transforms)
        np.testing.assert_allclose(ring.points_m, undeclared.points_m)

    def test_cli_aliases(self):
        args = plateau.build_parser().parse_args(['stl', '--root', 'input', '--out', 'out.stl',
            '--west', '139', '--east', '140', '--south', '35', '--north', '36',
            '--mm-per-m', '.1', '--export-origin', 'origin.json'])
        self.assertEqual(args.scale, .1)
        self.assertEqual(args.origin_out, 'origin.json')

    def test_stream_memory_and_separate_tiles_share_origin(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / 'input'; source.mkdir()
            # The CRS is inherited from the CityModel envelope, not the posList.
            for i in range(2):
                east, north = -6155 + i * 120, -34454
                triples = [(north, east, 1.8), (north, east + 12, 1.8),
                           (north + 12, east, 3.6)]
                content = ' '.join(' '.join(map(str, p)) for p in triples)
                (source / f'bldg_{i}.gml').write_text(
                    '<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0" '
                    'xmlns:gml="http://www.opengis.net/gml"><gml:boundedBy>'
                    '<gml:Envelope srsName="EPSG:6677" srsDimension="3"/>'
                    '</gml:boundedBy><gml:Polygon><gml:exterior><gml:LinearRing>'
                    f'<gml:posList>{content}</gml:posList></gml:LinearRing></gml:exterior>'
                    '</gml:Polygon></core:CityModel>')
            origin = ModelCoordinates(-6155, -34454, mm_per_m=1 / 120, vertical_exaggeration=200 / 3)
            origin.write(root / 'origin.json')
            args = argparse.Namespace(root=source, out=root / 'stream.stl', west=139, east=140,
                south=35, north=36, max_files=None, file_offset=0, skip_files_larger_than_mb=None,
                lod='any', scale=.001, base_thickness=0, repair=False,
                origin_json=root / 'origin.json', origin_out=root / 'used_origin.json',
                vertical_exaggeration=1, z_offset=0)
            with contextlib.redirect_stderr(io.StringIO()):
                plateau.build_stl(args)
                streamed = stl_triangles(args.out)
                args.out = root / 'memory.stl'
                plateau.build_stl_in_memory(args)
                np.testing.assert_array_equal(stl_triangles(args.out), streamed)
                parts = []
                for i in range(2):
                    args.file_offset, args.max_files = i, 1
                    args.out = root / f'part_{i}.stl'
                    plateau.build_stl(args)
                    parts.append(stl_triangles(args.out))
            np.testing.assert_array_equal(np.concatenate(parts), streamed)
            np.testing.assert_allclose(streamed[:, 0], [[0, 0, 1], [1, 0, 1]], atol=1e-7)
            self.assertEqual(ModelCoordinates.read(root / 'used_origin.json'), origin)

    def test_legacy_auto_center_is_retained(self):
        args = argparse.Namespace(scale=.1)
        bounds = plateau.Bounds(10, 20, 30, 20, 40, 50)
        c = plateau.coordinates_for_bounds(args, bounds)
        self.assertEqual(c.point(bounds.center()), (0, 0, 0))


if __name__ == '__main__':
    unittest.main()

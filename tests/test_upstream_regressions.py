"""Regressions for the actual upstream sweep and normalization entry points."""
import math
import unittest

import numpy as np

from subway_n02_welded_tube_model import make_welded_tube, make_hollow_tube, _sweep_frames
from yamanote_n02_overlay import make_tube
from subway_pipeline import ROUTES, make_alias_index, match_feature_route, merged_confidence


class SweepTests(unittest.TestCase):
    def assert_closed_positive(self, mesh):
        self.assertTrue(mesh.is_watertight)
        self.assertTrue(mesh.is_winding_consistent)
        self.assertGreater(mesh.volume, 0.)
        self.assertTrue(np.isfinite(mesh.vertices).all())

    def test_straight_solid_winding_and_polygon_volume(self):
        mesh = make_welded_tube([[0, 0, 0], [0, 0, 10]], 2., 24)
        self.assert_closed_positive(mesh)
        self.assertAlmostEqual(mesh.volume, .5 * 24 * math.sin(2 * math.pi / 24) * 4 * 10, places=10)

    def test_bent_solid_winding(self):
        mesh = make_welded_tube([[0, 0, 0], [1, 0, 10], [4, 2, 20]], .5, 16)
        self.assert_closed_positive(mesh)

    def test_repeated_endpoint_loop_has_indexed_seam(self):
        angles = np.linspace(0, 2 * np.pi, 41)
        points = np.c_[10 * np.cos(angles), 10 * np.sin(angles), angles * 0]
        mesh = make_welded_tube(points, .5, 16)
        self.assert_closed_positive(mesh)
        self.assertEqual(len(mesh.vertices), 40 * 16)
        self.assertEqual(mesh.euler_number, 0)

    def test_nonplanar_loop(self):
        a = np.linspace(0, 2 * np.pi, 61)
        p = np.c_[10 * np.cos(a), 10 * np.sin(a), np.sin(3 * a)]
        self.assert_closed_positive(make_welded_tube(p, .4, 16))

    def test_yamanote_uses_shared_closed_sweep(self):
        points = [[0, 0, 0], [20, 0, 0], [20, 20, 0], [0, 20, 0], [0, 0, 0]]
        actual = make_tube(points, .5, 12)
        expected = make_welded_tube(points, .5, 12)
        self.assert_closed_positive(actual)
        np.testing.assert_array_equal(actual.vertices, expected.vertices)
        np.testing.assert_array_equal(actual.faces, expected.faces)

    def test_open_hollow_has_annular_rims_and_no_bore_disc(self):
        mesh = make_hollow_tube([[0, 0, 0], [0, 0, 10]], 2., 1., 32)
        self.assert_closed_positive(mesh)
        expected = .5 * 32 * math.sin(2 * math.pi / 32) * (4 - 1) * 10
        self.assertAlmostEqual(mesh.volume, expected, places=10)
        # A ray along the bore hits no wall; one through the annulus hits both rims.
        self.assertFalse(mesh.ray.intersects_any([[0, 0, -1]], [[0, 0, 1]])[0])
        self.assertTrue(mesh.ray.intersects_any([[1.5, 0, -1]], [[0, 0, 1]])[0])

    def test_closed_hollow(self):
        a = np.linspace(0, 2 * np.pi, 51)
        p = np.c_[10 * np.cos(a), 10 * np.sin(a), a * 0]
        self.assert_closed_positive(make_hollow_tube(p, .8, .4, 16, closed=True))

    def test_duplicate_points_and_degenerate_route(self):
        mesh = make_welded_tube([[0, 0, 0], [0, 0, 0], [0, 0, 5]], 1., 8)
        self.assert_closed_positive(mesh)
        self.assertIsNone(make_welded_tube([[0, 0, 0], [0, 0, 0]], 1., 8))

    def test_frames_do_not_flip_near_vertical(self):
        z = np.linspace(0, 20, 81)
        p = np.c_[.25 * np.sin(z / 4), .15 * np.cos(z / 4), z]
        _, normal, _, _ = _sweep_frames(p, False)
        self.assertGreater(float(np.einsum('ij,ij->i', normal[1:], normal[:-1]).min()), .99)

    def test_invalid_parameters_rejected(self):
        p = [[0, 0, 0], [0, 0, 1]]
        for radius, sides in [(0, 12), (1, 2), (float('nan'), 12)]:
            with self.assertRaises(ValueError):
                make_welded_tube(p, radius, sides)
        with self.assertRaises(ValueError):
            make_hollow_tube(p, 1., 1., 12)


class RouteTests(unittest.TestCase):
    def setUp(self):
        self.aliases = make_alias_index(ROUTES)

    def test_shinjuku_is_not_hibiya(self):
        self.assertEqual(match_feature_route({'operator': 'Toei', 'line': 'Shinjuku Line'}, self.aliases),
                         ('Toei Subway', 'Shinjuku Line'))

    def test_route_codes_are_full_typed_field_matches(self):
        self.assertIsNone(match_feature_route({'description': 'building', 'line': 'unknown'}, self.aliases))
        self.assertIsNone(match_feature_route({'operator': 'Toei', 'line': 'G'}, self.aliases))
        self.assertEqual(match_feature_route({'operator': 'Tokyo Metro', 'route_short_name': 'G'}, self.aliases),
                         ('Tokyo Metro', 'Ginza Line'))

    def test_confidence_is_weakest_including_unknown(self):
        self.assertEqual(merged_confidence([{'confidence': 'official'}, {'confidence': 'estimated'}]), 'estimated')
        self.assertEqual(merged_confidence([{'confidence': 'official'}, {'confidence': 'future_unknown_value'}]), 'unknown')
        self.assertEqual(merged_confidence([]), 'unknown')


if __name__ == '__main__':
    unittest.main()

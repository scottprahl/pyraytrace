"""Pytest coverage for :mod:`pyraytrace.pyraytrace`."""

# ruff: noqa: D103

import numpy as np
import pytest

import pyraytrace.pyraytrace as geom


def test_plane_distance_and_parallel_case():
    """Plane distance should return finite hit distance or infinity if parallel."""
    plane = geom.Plane(xyz=(0, 0, 5), uvw=(0, 0, 1))
    ray_hit = geom.Ray(xyz=(0, 0, 0), uvw=(0, 0, 2))
    ray_parallel = geom.Ray(xyz=(0, 0, 0), uvw=(1, 0, 0))

    assert plane.distance(ray_hit) == pytest.approx(2.5)
    assert plane.distance(ray_parallel) == np.inf


def test_plane_is_in_plane_uses_tolerance():
    """Plane membership should honor the built-in tolerance threshold."""
    plane = geom.Plane(xyz=(0, 0, 1), uvw=(0, 0, 1))
    assert plane.is_in_plane(np.array([2, -3, 1 + 1e-7]))
    assert not plane.is_in_plane(np.array([0, 0, 1 + 1e-4]))


def test_plane_draw_zy_covers_both_equation_branches(monkeypatch):
    """draw_zy should handle both uvw[2] and uvw[1] branch equations."""
    calls = []
    monkeypatch.setattr(geom.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    plane1 = geom.Plane(xyz=(0, 0, 2), uvw=(0, 1, 1))
    plane1.draw_zy(ymin=0, ymax=2, zmin=0, zmax=5)

    plane2 = geom.Plane(xyz=(0, 1, 0), uvw=(0, 1, 0))
    plane2.draw_zy(ymin=0, ymax=3, zmin=0, zmax=2)

    assert len(calls) == 2
    assert np.allclose(calls[0][0][0], [2, 0])
    assert np.allclose(calls[0][0][1], [0, 2])
    assert calls[0][0][2] == "k"
    assert np.allclose(calls[1][0][0], [0, 2])
    assert np.allclose(calls[1][0][1], [1, 1])
    assert calls[1][0][2] == "k"


def test_ray_reflect_and_move_with_drawing(monkeypatch):
    """Ray reflection should update direction, and move should update position."""
    calls = []
    monkeypatch.setattr(geom.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    ray = geom.Ray(xyz=(0, 0, 0), uvw=(0, 1, -1))
    plane = geom.Plane(xyz=(0, 0, 0), uvw=(0, 0, 1))
    ray.reflect_from_plane(plane)
    assert np.allclose(ray.uvw, [0, 1, 1])

    ray.move(2, draw_zy=True)
    assert np.allclose(ray.xyz, [0, 2, 2])
    assert calls[0][0] == ([0, 2], [0, 2], "b")


def test_refract_handles_normal_flip_and_reflection_branch():
    """refract should support normal flipping and total internal reflection."""
    straight = geom.refract(
        uvw=np.array([0.0, 0.0, -1.0]),
        normal=np.array([0.0, 0.0, 1.0]),
        ni=1.0,
        nt=1.5,
    )
    assert np.allclose(straight, [0, 0, -1])

    reflected = geom.refract(
        uvw=np.array([np.sqrt(3) / 2, 0.0, 0.5]),
        normal=np.array([0.0, 0.0, 1.0]),
        ni=1.0,
        nt=1.5,
    )
    assert np.allclose(reflected, [np.sqrt(3) / 2, 0.0, -0.5])


def test_sphere_distance_covers_hit_inside_and_miss():
    """Sphere distance should cover external hit, internal hit, and miss cases."""
    sphere = geom.Sphere(xyz=(0, 0, 0), R=1, n=1.4)

    ray_hit = geom.Ray(xyz=(0, 0, -3), uvw=(0, 0, 1))
    assert sphere.distance(ray_hit) == pytest.approx(2.0)

    ray_inside = geom.Ray(xyz=(0, 0, 0), uvw=(0, 0, 1))
    assert sphere.distance(ray_inside) == pytest.approx(1.0)

    ray_miss = geom.Ray(xyz=(0, 0, -3), uvw=(0, 1, 0))
    assert sphere.distance(ray_miss) == np.inf


def test_sphere_distance_returns_inf_when_surface_intersection_is_behind():
    """Sphere distance should return inf when the intersection is behind the ray."""
    sphere = geom.Sphere(xyz=(0, 0, 0), R=1, n=1.4)
    ray_outward_from_surface = geom.Ray(xyz=(0, 0, 1), uvw=(0, 0, 1))
    assert sphere.distance(ray_outward_from_surface) == np.inf


def test_sphere_unit_normal_and_refract_delegation(monkeypatch):
    """Sphere refract should delegate to module refract with proper indices."""
    sphere = geom.Sphere(xyz=(1, 2, 3), R=2, n=1.6)
    point = np.array([3.0, 2.0, 3.0])
    normal = sphere.unit_normal_at(point)
    assert np.allclose(normal, [1, 0, 0])

    calls = []

    def fake_refract(uvw, normal_vec, ni, nt):
        """Capture refract call arguments while returning a sentinel vector."""
        calls.append((uvw.copy(), normal_vec.copy(), ni, nt))
        return np.array([9.0, 8.0, 7.0])

    monkeypatch.setattr(geom, "refract", fake_refract)
    ray = geom.Ray(xyz=point, uvw=(0, 0, 1))

    outside = sphere.refract(ray, outside=True)
    inside = sphere.refract(ray, outside=False)

    assert np.allclose(outside, [9, 8, 7])
    assert np.allclose(inside, [9, 8, 7])
    assert calls[0][2:] == (1, 1.6)
    assert calls[1][2:] == (1.6, 1)


def test_sphere_draw_zy_respects_side_selection(monkeypatch):
    """Sphere draw should emit traces according to side selection."""
    calls = []
    monkeypatch.setattr(geom.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    sphere = geom.Sphere(xyz=(0, 0, 0), R=2)
    sphere.draw_zy(side="both")
    sphere.draw_zy(side="right")
    sphere.draw_zy(side="left")

    assert len(calls) == 4


def test_prism_unit_normal_distance_draw_and_refract(monkeypatch):
    """Prism helpers should resolve normals, distance, drawing, and refraction."""
    A = geom.Plane(xyz=(0, 0, 1), uvw=(0, 0, 1))
    B = geom.Plane(xyz=(0, 1, 0), uvw=(0, 1, 0))
    C = geom.Plane(xyz=(0, -1, 0), uvw=(0, -1, 0))
    prism = geom.Prism(A, B, C, n=1.5)

    assert np.allclose(prism.unit_normal_at(np.array([0, 0, 1])), [0, 0, 1])
    assert np.allclose(prism.unit_normal_at(np.array([0, 1, 4])), [0, 1, 0])
    assert np.allclose(prism.unit_normal_at(np.array([0, -1, 4])), [0, -1, 0])
    assert np.allclose(prism.unit_normal_at(np.array([0, 0, 0])), [0, 0, 0])

    ray = geom.Ray(xyz=(0, 0, 0), uvw=(0, 0, 1))
    assert prism.distance(ray) == pytest.approx(1.0)

    side_calls = []
    monkeypatch.setattr(A, "draw_zy", lambda **kwargs: side_calls.append(("A", kwargs)))
    monkeypatch.setattr(B, "draw_zy", lambda **kwargs: side_calls.append(("B", kwargs)))
    monkeypatch.setattr(C, "draw_zy", lambda **kwargs: side_calls.append(("C", kwargs)))
    prism.draw_zy(ymin=-2, ymax=2, zmin=0, zmax=5)
    assert [name for name, _ in side_calls] == ["A", "B", "C"]

    refract_calls = []

    def fake_refract(uvw, normal_vec, ni, nt):
        """Capture refract arguments and return a sentinel direction."""
        refract_calls.append((uvw.copy(), normal_vec.copy(), ni, nt))
        return np.array([1.0, 0.0, 0.0])

    monkeypatch.setattr(geom, "refract", fake_refract)
    prism.refract(geom.Ray(xyz=(0, 0, 1), uvw=(0, 0, 1)), outside=True)
    prism.refract(geom.Ray(xyz=(0, 0, 1), uvw=(0, 0, 1)), outside=False)
    assert refract_calls[0][2:] == (1, 1.5)
    assert refract_calls[1][2:] == (1.5, 1)


def test_prism_distance_replaces_non_forward_intersections_with_large_value():
    """Prism distance should map non-forward intersections to the sentinel 999."""
    class _Face:
        """Minimal face stub with controllable intersection distance."""

        def __init__(self, d):
            """Store the synthetic distance value."""
            self._d = d
            self.uvw = np.array([0, 0, 1])

        def distance(self, _ray):
            """Return the configured distance for any ray."""
            return self._d

        def is_in_plane(self, _point):
            """Report no matches for synthetic in-plane checks."""
            return False

    prism = geom.Prism(_Face(0), _Face(-3), _Face(np.inf), n=1.4)
    assert prism.distance(geom.Ray()) == 999


def test_lens_routes_distance_refraction_and_drawing_calls():
    """Lens should route distance/refraction/draw calls to the right surface."""
    class _Surface:
        """Minimal lens surface stub used to verify delegation logic."""

        def __init__(self, d):
            """Store distance and call tracking containers."""
            self.d = d
            self.refract_calls = []
            self.draw_calls = []

        def distance(self, _ray):
            """Return configured distance for the provided ray."""
            return self.d

        def refract(self, _ray, value):
            """Track refract factor and return a synthetic vector."""
            self.refract_calls.append(value)
            return np.array([value, 0, 0])

        def draw_zy(self, side=None):
            """Record which side draw call was requested."""
            self.draw_calls.append(side)

    s1 = _Surface(2.0)
    s2 = _Surface(5.0)
    lens = geom.Lens(s1, s2, refractive_index=1.5, thickness=3.0)
    ray = geom.Ray()

    assert lens.distance(ray, which_surface=1) == pytest.approx(2.0)
    assert lens.distance(ray, which_surface=2) == pytest.approx(5.0)
    assert np.allclose(lens.refract(ray, which_surface=1), [2 / 3, 0, 0])
    assert np.allclose(lens.refract(ray, which_surface=2), [1.5, 0, 0])
    assert s1.refract_calls == [2 / 3]
    assert s2.refract_calls == [1.5]

    lens.draw_zy()
    assert s1.draw_calls == ["left"]
    assert s2.draw_calls == ["right"]


def test_thin_lens_distance_refraction_and_drawing(monkeypatch):
    """Thin lens helpers should compute distance, scalar refraction, and draw."""
    calls = []
    monkeypatch.setattr(geom.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    lens = geom.ThinLens(focal_length=10, vertex=5, diameter=8)
    ray = geom.Ray(xyz=(0, 2, 1), uvw=(0, 0, 2))
    assert lens.distance(ray) == pytest.approx(2.0)

    parallel = geom.Ray(xyz=(0, 2, 1), uvw=(1, 0, 0))
    assert lens.distance(parallel) == np.inf

    expected = 10 / np.sqrt((10 / 2 - 2) ** 2 + 10**2)
    assert lens.refract(ray) == pytest.approx(expected)

    lens.draw_zy()
    assert calls[0][0] == ([5, 5], [-4.0, 4.0], ":r")

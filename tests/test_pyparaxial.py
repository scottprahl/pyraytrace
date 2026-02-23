"""Pytest coverage for :mod:`pyraytrace.pyparaxial`."""

# ruff: noqa: D103

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg", force=True)

import pyraytrace.pyparaxial as paraxial


def test_ray_aim_methods_update_height_and_slope():
    ray = paraxial.Ray()
    ray.aim(z0=1, y0=2, z1=3, y1=6)
    assert np.allclose(ray.ynu, [2, 2])

    ray.aim2((0, -1), (2, 3))
    assert np.allclose(ray.ynu, [-1, 2])


def test_ray_string_representations_include_values():
    ray = paraxial.Ray(y=1.25, nu=-0.5)
    assert "Ray y=" in str(ray)
    assert "nu=" in str(ray)
    assert "Ray(y=" in repr(ray)


def test_paraxial_element_draw_uses_vertical_segment(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    element = paraxial.ParaxialElement(z=2, h=3, color="green")
    element.draw()

    assert calls == [(([2, 2], [-3, 3]), {"color": "green"})]


def test_mirror_sets_expected_abcd_matrix_and_draws(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    mirror = paraxial.Mirror(z=1, h=2, R=10, color="cyan")
    assert np.allclose(mirror.ABCD, [[1, 0], [0.2, 1]])

    mirror.draw()
    assert calls == [(([1, 1], [-2, 2]), {"color": "cyan"})]


@pytest.mark.parametrize(
    ("focal_length", "expected_top", "expected_bottom"),
    [
        (5, [1.5, 2, 1.5], [-1.5, -2, -1.5]),
        (-5, [2.5, 2, 2.5], [-2.5, -2, -2.5]),
    ],
)
def test_thin_lens_draw_handles_converging_and_diverging_branches(
    monkeypatch, focal_length, expected_top, expected_bottom
):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    lens = paraxial.ThinLens(z=10, h=2, f=focal_length, width=0.25, color="blue")
    lens.draw()

    assert len(calls) == 3
    assert calls[0] == (([10, 10], [-2, 2]), {"color": "blue"})
    assert np.allclose(calls[1][0][1], expected_top)
    assert np.allclose(calls[2][0][1], expected_bottom)
    assert all(call[1]["color"] == "blue" for call in calls)


def test_plane_draw_uses_dotted_style(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    plane = paraxial.Plane(z=4, h=1.5, color="black")
    plane.draw()

    assert calls[0][0] == ([4, 4], [-1.5, 1.5], ":")
    assert calls[0][1] == {"color": "black"}


def test_aperture_draws_expected_four_segments(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    aperture = paraxial.Aperture(z=1, h=4, D=2, width=0.5, color="orange")
    aperture.draw()

    assert len(calls) == 4
    assert calls[0] == (([1, 1], [1, 4]), {"color": "orange"})
    assert calls[1] == (([0.0, 2.0], [1, 1]), {"color": "orange"})
    assert calls[2] == (([1, 1], [-1, -4]), {"color": "orange"})
    assert calls[3] == (([0.0, 2.0], [-1, -1]), {"color": "orange"})


def test_object_draw_uses_annotate_arrow(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "annotate", lambda *args, **kwargs: calls.append((args, kwargs)))

    obj = paraxial.Object(z=3, h=2, color="red")
    obj.draw()

    assert len(calls) == 1
    args, kwargs = calls[0]
    assert args == ("",)
    assert kwargs["xy"] == (3, 0)
    assert kwargs["xytext"] == (3, 2)
    assert kwargs["arrowprops"] == {"arrowstyle": "<-"}


def test_paraxial_system_minmax_and_text_representations():
    system = paraxial.ParaxialSystem()
    system.append(paraxial.Plane(z=5, h=1))
    system.append(paraxial.Mirror(z=-2, h=2, R=4))

    assert system.minmax() == (-2, 5)
    assert "Plane z=" in str(system)
    assert "Mirror(z=" in repr(system)


def test_draw_ray_through_system_plots_each_free_space_segment(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    system = paraxial.ParaxialSystem()
    system.append(paraxial.Plane(z=0, h=1))
    system.append(paraxial.ThinLens(z=2, h=1, f=4))
    system.append(paraxial.Mirror(z=5, h=1, R=10))

    ray = paraxial.Ray(y=1, nu=0.5)
    system.draw_ray_through_system(ray, color="magenta")

    assert len(calls) == 2
    assert calls[0] == (([0, 2], [1.0, 2.0]), {"color": "magenta"})
    assert calls[1] == (([2, 5], [2.0, 2.0]), {"color": "magenta"})


def test_paraxial_system_draw_plots_axis_and_calls_element_draw(monkeypatch):
    calls = []
    monkeypatch.setattr(paraxial.plt, "plot", lambda *args, **kwargs: calls.append((args, kwargs)))

    class _Element:
        def __init__(self, z):
            self.z = z
            self.draw_calls = 0

        def draw(self):
            self.draw_calls += 1

    e1 = _Element(-1)
    e2 = _Element(3)
    system = paraxial.ParaxialSystem()
    system.append(e1)
    system.append(e2)

    system.draw()

    assert calls[0][0] == ([-1, 3], [0, 0], "k")
    assert e1.draw_calls == 1
    assert e2.draw_calls == 1

# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
"""Paraxial (2D) ray-tracing primitives.

This module provides simple paraxial elements and ray propagation helpers for
visualizing optical systems in the ``z-y`` plane.
"""

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "Ray",
    "ParaxialElement",
    "Mirror",
    "ThinLens",
    "Plane",
    "Aperture",
    "Object",
    "ParaxialSystem",
]


class Ray:
    """A 2D ray represented by height and slope."""

    def __init__(self, y=0, nu=1):
        """Initialize a paraxial ray.

        Args:
            y: Initial transverse height.
            nu: Initial slope angle in radians.
        """
        self.ynu = np.array([y, nu], dtype=float)

    def __str__(self):
        """Return a compact human-readable representation."""
        return "Ray y=%.3f nu=%.3f" % (self.ynu[0], self.ynu[1])

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "Ray(y=%f, nu=%f)" % (self.ynu[0], self.ynu[1])

    def aim(self, z0=0, y0=0, z1=1, y1=1):
        """Set the ray to pass through two points.

        Args:
            z0: Axial coordinate of the first point.
            y0: Height coordinate of the first point.
            z1: Axial coordinate of the second point.
            y1: Height coordinate of the second point.
        """
        self.ynu[0] = y0
        self.ynu[1] = (y1 - y0) / (z1 - z0)

    def aim2(self, p, q):
        """Set the ray to pass through two ``(z, y)`` points.

        Args:
            p: First ``(z, y)`` point.
            q: Second ``(z, y)`` point.
        """
        self.ynu[0] = p[1]
        self.ynu[1] = (q[1] - p[1]) / (q[0] - p[0])


class ParaxialElement:
    """A 2D optical element represented by an ABCD matrix."""

    def __init__(self, z=0, h=1, color="black"):
        """Initialize a paraxial element.

        Args:
            z: Axial position.
            h: Half-height used when drawing the element.
            color: Matplotlib color used when drawing.
        """
        self.z = z
        self.h = h
        self.ABCD = np.array([[1, 0], [0, 1]], dtype=float)
        self.color = color

    def __str__(self):
        """Return a compact human-readable representation."""
        return "ParaxialElement z=%.3f, h=%.3f" % (self.z, self.h)

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "ParaxialElement(z=%f, h=%f)" % (self.z, self.h)

    def draw(self):
        """Draw the element in the ``z-y`` plane."""
        plt.plot([self.z, self.z], [-self.h, self.h], color=self.color)


class Mirror(ParaxialElement):
    """A paraxial spherical mirror."""

    def __init__(self, z=0, h=1, R=np.inf, color="purple"):
        """Initialize a mirror.

        Args:
            z: Axial position of the mirror.
            h: Half-height used when drawing.
            R: Radius of curvature.
            color: Matplotlib color used when drawing.
        """
        super().__init__(z=z, h=h, color=color)
        self.R = R
        self.ABCD = np.array([[1, 0], [2 / R, 1]], dtype=float)

    def __str__(self):
        """Return a compact human-readable representation."""
        return "Mirror z=%.3f h=%.3f R=%.3f" % (self.z, self.h, self.R)

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "Mirror(z=%f, h=%f, R=%f)" % (self.z, self.h, self.R)

    def draw(self):
        """Draw the mirror in the ``z-y`` plane."""
        plt.plot([self.z, self.z], [-self.h, self.h], color=self.color)


class ThinLens(ParaxialElement):
    """A paraxial thin lens."""

    def __init__(self, z=0, h=1, f=10, width=0.2, color="blue"):
        """Initialize a thin lens.

        Args:
            z: Axial position of the lens.
            h: Half-height used when drawing.
            f: Focal length.
            width: Fractional width used for the drawn arrowhead.
            color: Matplotlib color used when drawing.
        """
        super().__init__(z=z, h=h, color=color)
        self.f = f
        self.width = width
        self.ABCD = np.array([[1, 0], [-1 / f, 1]], dtype=float)

    def __str__(self):
        """Return a compact human-readable representation."""
        return "ThinLens z=%.3f h=%.3f f=%.3f" % (self.z, self.h, self.f)

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "ThinLens(z=%f, h=%f, f=%f)" % (self.z, self.h, self.f)

    def draw(self):
        """Draw the thin lens in the ``z-y`` plane."""
        w = self.h * self.width
        plt.plot([self.z, self.z], [-self.h, self.h], color=self.color)
        if self.f > 0:
            plt.plot(
                [self.z - w, self.z, self.z + w],
                [self.h - w, self.h, self.h - w],
                color=self.color,
            )
            plt.plot(
                [self.z - w, self.z, self.z + w],
                [-self.h + w, -self.h, -self.h + w],
                color=self.color,
            )
        else:
            plt.plot(
                [self.z - w, self.z, self.z + w],
                [self.h + w, self.h, self.h + w],
                color=self.color,
            )
            plt.plot(
                [self.z - w, self.z, self.z + w],
                [-self.h - w, -self.h, -self.h - w],
                color=self.color,
            )


class Plane(ParaxialElement):
    """A paraxial reference plane."""

    def __init__(self, z=0, h=1, color="black"):
        """Initialize a plane element.

        Args:
            z: Axial position of the plane.
            h: Half-height used when drawing.
            color: Matplotlib color used when drawing.
        """
        super().__init__(z=z, h=h, color=color)

    def __str__(self):
        """Return a compact human-readable representation."""
        return "Plane z=%.3f h=%.3f" % (self.z, self.h)

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "Plane(z=%f, h=%f)" % (self.z, self.h)

    def draw(self):
        """Draw the plane in the ``z-y`` plane."""
        plt.plot([self.z, self.z], [-self.h, self.h], ":", color=self.color)


class Aperture(Plane):
    """A paraxial aperture stop."""

    def __init__(self, z=0, h=1, D=1, width=0.3, color="black"):
        """Initialize an aperture stop.

        Args:
            z: Axial position of the aperture.
            h: Total half-height used for drawing.
            D: Aperture diameter.
            width: Fractional bar width used when drawing the stop.
            color: Matplotlib color used when drawing.
        """
        super().__init__(z=z, h=h, color=color)
        self.D = D
        self.width = width

    def __str__(self):
        """Return a compact human-readable representation."""
        return "Aperture z=%.3f h=%.3f D=%.3f width=%.3f" % (
            self.z,
            self.h,
            self.D,
            self.width,
        )

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "Aperture(z=%f, h=%f, D=%f)" % (self.z, self.h, self.D)

    def draw(self):
        """Draw the aperture stop in the ``z-y`` plane."""
        r = self.D / 2
        w = self.width * self.h / 2
        plt.plot([self.z, self.z], [r, self.h], color=self.color)
        plt.plot([self.z - w, self.z + w], [r, r], color=self.color)
        plt.plot([self.z, self.z], [-r, -self.h], color=self.color)
        plt.plot([self.z - w, self.z + w], [-r, -r], color=self.color)


class Object(Plane):
    """A paraxial object marker."""

    def __init__(self, z=0, h=1, color="red"):
        """Initialize an object marker.

        Args:
            z: Axial position of the object.
            h: Object height.
            color: Matplotlib color used when drawing.
        """
        super().__init__(z=z, h=h, color=color)

    def __str__(self):
        """Return a compact human-readable representation."""
        return "Object z=%.3f h=%.3f" % (self.z, self.h)

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "Object(z=%f, h=%f)" % (self.z, self.h)

    def draw(self):
        """Draw the object as an arrow in the ``z-y`` plane."""
        plt.annotate(
            "",
            xy=(self.z, 0),
            xytext=(self.z, self.h),
            arrowprops={"arrowstyle": "<-"},
        )


class ParaxialSystem:
    """A sequence of paraxial elements."""

    def __init__(self):
        """Initialize an empty paraxial system."""
        self.elements = []

    def __str__(self):
        """Return a compact human-readable representation."""
        s = ""
        for e in self.elements:
            s += e.__str__() + "\n"
        return s

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        r = ""
        for e in self.elements:
            r += e.__repr__() + "\n"
        return r

    def append(self, elem):
        """Append an optical element to the system.

        Args:
            elem: Element to append.
        """
        self.elements.append(elem)

    def minmax(self):
        """Return the minimum and maximum element z locations.

        Returns:
            Tuple ``(zmin, zmax)`` across all elements.
        """
        zmin = 1e6
        zmax = -1e6
        for e in self.elements:
            zmin = min(zmin, e.z)
            zmax = max(zmax, e.z)
        return zmin, zmax

    def draw_ray_through_system(self, ray, color="black"):
        """Propagate and draw a paraxial ray through the system.

        Args:
            ray: Input paraxial ray.
            color: Matplotlib color for the plotted trajectory.
        """
        z_last = None
        ynu = np.copy(ray.ynu)

        for e in self.elements:
            if z_last is not None:
                prop = np.array([[1, e.z - z_last], [0, 1]])

                ynu_next = np.matmul(prop, ynu)
                plt.plot([z_last, e.z], [ynu[0], ynu_next[0]], color=color)

                ynu_bent = np.matmul(e.ABCD, ynu_next)
                ynu = np.copy(ynu_bent)
            z_last = e.z

    def draw(self):
        """Draw all elements in the system."""
        zmin, zmax = self.minmax()

        plt.plot([zmin, zmax], [0, 0], "k")
        for e in self.elements:
            e.draw()

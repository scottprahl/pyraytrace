# pylint: disable=invalid-name
"""Geometric (3D) ray-tracing primitives.

This module provides basic objects and helper routines for tracing rays through
planes, spheres, prisms, and lenses.
"""

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["Plane", "Ray", "Sphere", "Prism", "Lens", "ThinLens"]


class Plane:
    """A plane used as an optical surface."""

    def __init__(self, xyz=(0, 0, 0), uvw=(0, 0, 1)):
        """Initialize a plane.

        Args:
            xyz: A point on the plane in Cartesian coordinates.
            uvw: A unit normal vector for the plane.
        """
        self.xyz = np.array(xyz)
        self.uvw = np.array(uvw)
        self.D = -np.dot(xyz, uvw)

    def __str__(self):
        """Return a compact human-readable representation."""
        a = "xyz=[%.3f,%.3f,%.3f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "uvw=[%.3f,%.3f,%.3f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        length = np.dot(self.uvw, self.uvw)
        return a + ", " + b + " norm=%.4f" % length + " D=%f" % self.D

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "[%f,%f,%f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        return "Plane(" + a + ", " + b + ")"

    def draw_zy(self, ymin=0, ymax=1, zmin=0, zmax=1):
        """Draw the plane segment in the ``z-y`` view.

        Args:
            ymin: Minimum y value of the view window.
            ymax: Maximum y value of the view window.
            zmin: Minimum z value of the view window.
            zmax: Maximum z value of the view window.
        """
        if self.uvw[2] != 0:
            ymn = ymin
            ymx = ymax
            zmn = -(self.D + self.uvw[1] * ymin) / self.uvw[2]
            zmx = -(self.D + self.uvw[1] * ymax) / self.uvw[2]
            zmn = max(zmin, zmn)
            zmx = min(zmax, zmx)
            print("   zy=(%.2f,%.2f), zy=(%.2f,%.2f)" % (zmn, ymn, zmx, ymx))
            plt.plot([zmn, zmx], [ymn, ymx], "k")
            return

        if self.uvw[1] != 0:
            ymn = -(self.D + self.uvw[2] * zmin) / self.uvw[1]
            ymx = -(self.D + self.uvw[2] * zmax) / self.uvw[1]
            ymn = max(ymn, ymin)
            ymx = min(ymx, ymax)
            zmn = zmin
            zmx = zmax
            print("   zy=(%.2f,%.2f), zy=(%.2f,%.2f)" % (zmn, ymn, zmx, ymx))
            plt.plot([zmn, zmx], [ymn, ymx], "k")

    def distance(self, ray):
        """Return the signed intersection distance from a ray to the plane.

        Args:
            ray: Ray to test against the plane.

        Returns:
            Distance along the ray direction, or ``np.inf`` if parallel.
        """
        cos_angle = np.dot(ray.uvw, self.uvw)
        if abs(cos_angle) < 1e-8:
            return np.inf

        return -(np.dot(ray.xyz, self.uvw) + self.D) / cos_angle

    def is_in_plane(self, point):
        """Return whether a point lies on the plane.

        Args:
            point: Cartesian point to test.

        Returns:
            ``True`` when the point lies on the plane within tolerance.
        """
        dist = abs(np.dot(point, self.uvw) + self.D)
        return dist < 1e-6


class Ray:
    """A 3D ray represented by origin and direction cosines."""

    def __init__(self, xyz=(0, 0, 0), uvw=(0, 0, 1)):
        """Initialize a 3D ray.

        Args:
            xyz: Starting point of the ray.
            uvw: Direction vector of the ray.
        """
        self.xyz = np.array(xyz)
        self.uvw = np.array(uvw)

    def __str__(self):
        """Return a compact human-readable representation."""
        a = "xyz=[%.3f,%.3f,%.3f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "uvw=[%.3f,%.3f,%.3f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        length = np.dot(self.uvw, self.uvw)
        return a + ", " + b + " norm=%.4f" % length

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "[%f,%f,%f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        return "Ray(" + a + ", " + b + ")"

    def reflect_from_plane(self, plane):
        """Reflect the ray about a plane normal.

        Args:
            plane: Plane used as the reflecting surface.
        """
        a = np.dot(self.uvw, plane.uvw)
        out = self.uvw - 2 * a * plane.uvw
        self.uvw = out

    def move(self, d, draw_zy=False):
        """Advance the ray by a distance along its direction.

        Args:
            d: Distance to propagate.
            draw_zy: Whether to draw the segment in the ``z-y`` view.
        """
        dest = self.xyz + d * self.uvw
        if draw_zy:  # vertical is y and horizontal is z
            plt.plot([self.xyz[2], dest[2]], [self.xyz[1], dest[1]], "b")
        self.xyz = dest


def refract(uvw, normal, ni, nt):
    """Refract a direction vector using Snell's law.

    Args:
        uvw: Incoming unit direction vector.
        normal: Surface normal at the interaction point.
        ni: Incident refractive index.
        nt: Transmitted refractive index.

    Returns:
        Outgoing direction vector. On total internal reflection, returns the
        reflected direction.
    """
    cosine = np.dot(normal, uvw)
    if cosine < 0:
        cosine *= -1
        normal *= -1

    refractive_index = nt / ni
    a = refractive_index * cosine
    b = refractive_index**2 - 1
    disc = a**2 - b
    if disc < 0:  # reflected
        out = uvw - 2 * cosine * normal
    else:
        g = -a + np.sqrt(disc)
        out = refractive_index * uvw + g * normal

    return out


class Sphere:
    """A spherical optical surface."""

    def __init__(self, xyz=(0, 0, 0), R=1.0, n=1.0):
        """Initialize a sphere.

        Args:
            xyz: Center of the sphere in Cartesian coordinates.
            R: Sphere radius.
            n: Refractive index of the sphere medium.
        """
        self.xyz = np.array(xyz)
        self.radius = R
        self.n = n

    def __str__(self):
        """Return a compact human-readable representation."""
        a = "center=[%.3f,%.3f,%.3f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = ", radius = %f" % self.radius
        return a + b

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        return "Sphere(" + a + ", %f" % self.radius + ")"

    def draw_zy(self, ymax=np.inf, side="both"):
        """Draw the sphere cross-section in the ``z-y`` plane.

        Args:
            ymax: Maximum y value to draw.
            side: Which side to draw: ``"left"``, ``"right"``, or ``"both"``.
        """
        RR = np.sqrt(self.radius**2 - self.xyz[0] ** 2)
        yy = min(ymax, RR)
        y = np.linspace(-yy, yy, 50)
        r = RR**2 - (y - self.xyz[1]) ** 2
        np.place(r, r < 0, 0)
        z = np.sqrt(r)
        if side in ("both", "right"):
            plt.plot(z + self.xyz[2], y, "k")
        if side in ("both", "left"):
            plt.plot(-z + self.xyz[2], y, "k")

    def unit_normal_at(self, point):
        """Return the outward unit normal at a surface point.

        Args:
            point: Point on the sphere surface.

        Returns:
            Outward unit normal vector.
        """
        diff = point - self.xyz
        mag = np.sqrt(np.dot(diff, diff))
        return diff / mag

    def distance(self, ray):
        """Return the nearest positive intersection distance with the sphere.

        Args:
            ray: Ray to intersect with the sphere.

        Returns:
            Nearest positive distance, or ``np.inf`` when there is no
            forward intersection.
        """
        OS = ray.xyz - self.xyz
        b = 2 * np.dot(ray.uvw, OS)
        c = np.dot(OS, OS) - self.radius * self.radius

        disc = b * b - 4 * c
        if disc < 0:
            return np.inf

        disc = np.sqrt(disc)
        d1 = (-b - disc) / 2
        d2 = (-b + disc) / 2

        if d1 > 1e-6:
            return d1
        if d2 > 1e-6:
            return d2
        return np.inf

    def refract(self, ray, outside=True):
        """Compute the refracted ray direction at the sphere.

        Args:
            ray: Incoming ray at the sphere surface.
            outside: ``True`` when entering the sphere, ``False`` when exiting.

        Returns:
            Refracted (or reflected) direction vector.
        """
        normal = self.unit_normal_at(ray.xyz)
        if outside:
            return refract(ray.uvw, normal, 1, self.n)
        return refract(ray.uvw, normal, self.n, 1)


class Prism:
    """A prism defined by three planar faces."""

    def __init__(self, A, B, C, n):
        """Initialize a prism.

        Args:
            A: Plane for the first side.
            B: Plane for the second side.
            C: Plane for the third side.
            n: Refractive index inside the prism.
        """
        self.A = A
        self.B = B
        self.C = C
        self.n = n

    def __str__(self):
        """Return a compact human-readable representation."""
        return self.A.__str__() + "\n" + self.B.__str__() + "\n" + self.C.__str__()

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return self.A.__repl__() + self.B.__repl__() + self.C.__repl__()

    def draw_zy(self, ymin=0, ymax=1, zmin=0, zmax=1):
        """Draw prism faces in the ``z-y`` view.

        Args:
            ymin: Minimum y value of the view window.
            ymax: Maximum y value of the view window.
            zmin: Minimum z value of the view window.
            zmax: Maximum z value of the view window.
        """
        self.A.draw_zy(ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
        self.B.draw_zy(ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
        self.C.draw_zy(ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)

    def unit_normal_at(self, point):
        """Return the outward normal for the prism face at a point.

        Args:
            point: Candidate point on a prism face.

        Returns:
            Unit normal of the matching face, or a zero vector if unmatched.
        """
        if self.A.is_in_plane(point):
            return self.A.uvw
        if self.B.is_in_plane(point):
            return self.B.uvw
        if self.C.is_in_plane(point):
            return self.C.uvw

        # need to fail here
        return np.array([0, 0, 0])

    def distance(self, ray):
        """Return the nearest positive distance from a ray to a prism face.

        Args:
            ray: Ray to intersect with prism faces.

        Returns:
            Nearest positive face distance.
        """
        d1 = self.A.distance(ray)
        d2 = self.B.distance(ray)
        d3 = self.C.distance(ray)
        dd = np.array([d1, d2, d3])
        np.place(dd, dd <= 1e-6, 999)
        return min(dd)

    def refract(self, ray, outside):
        """Compute the refracted ray direction at a prism face.

        Args:
            ray: Incoming ray at a prism face.
            outside: ``True`` when entering, ``False`` when exiting.

        Returns:
            Refracted (or reflected) direction vector.
        """
        normal = self.unit_normal_at(ray.xyz)
        if outside:
            return refract(ray.uvw, normal, 1, self.n)
        return refract(ray.uvw, normal, self.n, 1)


class Lens:
    """A thick lens represented by two surfaces."""

    def __init__(self, surface1, surface2, refractive_index, thickness):
        """Initialize a lens.

        Args:
            surface1: First lens surface.
            surface2: Second lens surface.
            refractive_index: Lens refractive index.
            thickness: Physical lens thickness.
        """
        self.surface1 = surface1
        self.surface2 = surface2
        self.refractive_index = refractive_index
        self.thickness = thickness

    def __str__(self):
        """Return a compact human-readable representation."""
        a = str(self.surface1)
        b = str(self.surface2)
        c = "refractive index = %f" % self.refractive_index
        d = "thickness = %f" % self.thickness
        return a + "\n" + b + "\n" + c + "\n" + d

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        a = repr(self.surface1)
        b = repr(self.surface2)
        c = ", %f, %f)" % (self.refractive_index, self.thickness)
        return "Lens(" + a + "," + b + c

    def distance(self, ray, which_surface):
        """Return distance from a ray to one lens surface.

        Args:
            ray: Ray to intersect with lens surfaces.
            which_surface: Surface selector (`1` for first, otherwise second).

        Returns:
            Intersection distance to the selected surface.
        """
        if which_surface == 1:
            return self.surface1.distance(ray)

        return self.surface2.distance(ray)

    def refract(self, ray, which_surface):
        """Refract a ray at one lens surface.

        Args:
            ray: Incoming ray at the chosen surface.
            which_surface: Surface selector (`1` for first, otherwise second).

        Returns:
            Refracted (or reflected) direction vector.
        """
        if which_surface == 1:
            return self.surface1.refract(ray, 1 / self.refractive_index)

        return self.surface2.refract(ray, self.refractive_index)

    def draw_zy(self):
        """Draw the lens surfaces in the ``z-y`` plane."""
        self.surface1.draw_zy(side="left")
        self.surface2.draw_zy(side="right")


class ThinLens:
    """A thin-lens approximation."""

    def __init__(self, focal_length, vertex, diameter=10):
        """Initialize a thin lens.

        Args:
            focal_length: Lens focal length.
            vertex: Axial location of the lens plane.
            diameter: Lens diameter used for drawing.
        """
        self.f = focal_length
        self.vertex = vertex
        self.diameter = diameter

    def __str__(self):
        """Return a compact human-readable representation."""
        a = "focal length = %f" % self.f
        b = "vertex = %f" % self.vertex
        c = "diameter = %f" % self.diameter
        return a + "\n" + b + "\n" + c

    def __repr__(self):
        """Return an unambiguous representation for debugging."""
        return "ThinLens(%f, %f, diameter=%f)" % (self.f, self.vertex, self.diameter)

    def distance(self, ray):
        """Return distance from a ray origin to the thin lens plane.

        Args:
            ray: Ray to intersect with the thin lens plane.

        Returns:
            Distance to the lens plane, or ``np.inf`` when parallel.
        """
        z0 = ray.xyz[2]
        z1 = self.vertex
        w = ray.uvw[2]

        if w == 0:
            return np.inf
        return (z1 - z0) / w

    def refract(self, ray):
        """Compute a simplified refraction factor for the thin lens model.

        Args:
            ray: Incoming ray.

        Returns:
            Scalar factor used by the legacy thin-lens model.
        """
        w = ray.uvw[2]
        yf = self.f / w
        y0 = ray.xyz[1]
        return self.f / np.sqrt((yf - y0) ** 2 + self.f**2)

    def draw_zy(self):
        """Draw the thin lens in the ``z-y`` plane."""
        plt.plot([self.vertex, self.vertex], [-self.diameter / 2, self.diameter / 2], ":r")

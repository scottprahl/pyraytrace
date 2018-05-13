# pylint: disable=invalid-name
"""
Basic routines for ray tracing

To do:
    * Add lenses and mirrors
    * Properly document

Scott Prahl
May 2018
"""

import numpy as np
import matplotlib.pyplot as plt

__all__ = ['Plane',
           'Ray',
           'Sphere']

class Plane(object):
    """
    A class to help to ray-trace planar objects

    A plane is defined by the equation:

        u*x + v*y + w*z + D = 0

    where (u,v,w) is a unit normal vector to the plane.
    """

    def __init__(self, xyz=(0, 0, 0), uvw=(0, 0, 1)):
        """
        Args:
            xyz: array describing any point in the plane
            uvw: array with the direction cosines for the unit vector normal to the plane
        """

        self.xyz = np.array(xyz)
        self.uvw = np.array(uvw)
        self.D = -np.dot(xyz, uvw)

    def __str__(self):
        a = "xyz=[%.2f,%.2f,%.2f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "uvw=[%.2f,%.2f,%.2f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        length = np.dot(self.uvw, self.uvw)
        return a + ", " + b + " norm=%.4f" % length + " D=%f" % self.D

    def __repr__(self):
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "[%f,%f,%f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        return "Plane(" + a + ", " + b + ")"

    def draw_zy(self, ymin=0, ymax=1):
        """
        Draw representation in the zy-plane
        """
        zmin = -(self.D + self.uvw[1] * ymin) / self.uvw[2]
        zmax = -(self.D + self.uvw[1] * ymax) / self.uvw[2]
        plt.plot([zmin, zmax], [ymin, ymax], 'k')

    def distance(self, ray):
        """
        distance from start of ray to sphere
        """
        cos_angle = np.dot(ray.uvw, self.uvw)
        if abs(cos_angle) < 1e-8:
            return np.inf
        else:
            return -(np.dot(ray.xyz, self.uvw) + self.D) / cos_angle


class Ray(object):
    """
    A 3D ray specified by a starting point and a set of direction cosines
    """

    def __init__(self, xyz=(0, 0, 0), uvw=(0, 0, 1)):
        """
        Args:
            xyz: array describing the starting point for the ray
            uvw: array with the direction cosines
        """
        self.xyz = np.array(xyz)
        self.uvw = np.array(uvw)

    def __str__(self):
        a = "xyz=[%.2f,%.2f,%.2f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "uvw=[%.2f,%.2f,%.2f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        length = np.dot(self.uvw, self.uvw)
        return a + ", " + b + " norm=%.4f" % length

    def __repr__(self):
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "[%f,%f,%f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        return "Ray(" + a + ", " + b + ")"

    def reflect_from_plane(self, plane):
        """
        Spencer and Murty equation 46
        """
        a = np.dot(self.uvw, plane.uvw)
        out = self.uvw - 2 * a * plane.uvw
        self.uvw = out

    def move(self, d, draw_zy=False):
        """
        Spencer and Murty equation 5
        """
        dest = self.xyz + d * self.uvw
        if draw_zy:  # vertical is y and horizontal is z
            plt.plot([self.xyz[2], dest[2]], [self.xyz[1], dest[1]], 'b')
        self.xyz = dest


class Sphere(object):
    """
    A class to help to ray-trace spherical objects

    A sphere is defined by the equation:

        (x-x0)**2 + (y-y0)**2 + (z-z0)**2 = R**2

    where (x0,y0,z0) is the center of the sphere and R is the radius
    """

    def __init__(self, xyz=(0, 0, 0), R=1.0):
        """
        Args:
            xyz: array describing the center of the sphere
            R: radius of the sphere
        """

        self.xyz = np.array(xyz)
        self.radius = R

    def __str__(self):
        a = "center=[%.2f,%.2f,%.2f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = ", radius = %f" % self.radius
        return a + b

    def __repr__(self):
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        return "Sphere(" + a + ", %f" % self.radius + ")"

    def draw_zy(self):
        """
        Draw representation in the zy-plane
        """
        RR = np.sqrt(self.radius**2 - self.xyz[0]**2)
        y = np.linspace(-RR, RR, 50)
        z = np.sqrt(RR**2 - (y - self.xyz[1])**2)
        plt.plot(z + self.xyz[2], y, 'k')
        plt.plot(-z + self.xyz[2], y, 'k')

    def normal(self, point):
        """
        Return outward normal to point on sphere
        """
        diff = point - self.xyz
        mag = np.sqrt(np.dot(diff, diff))
        return diff / mag

    def distance(self, ray):
        """
        Return the distance of a ray to the sphere
        """
        OS = ray.xyz - self.xyz
        b = 2 * np.dot(ray.uvw, OS)
        c = np.dot(OS, OS) - self.radius * self.radius

        disc = b * b - 4 * c
        if disc < 0:
            return np.inf
        disc = np.sqrt(disc)
        d = (-b - disc) / 2
        if d > 1e-6:
            return d
        d = (-b + disc) / 2
        if d > 1e-6:
            return d
        return np.inf

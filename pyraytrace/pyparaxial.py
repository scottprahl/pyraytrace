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
           'Sphere',
           'Prism',
           'Lens',
           'ThinLens']


class Ray:
    """
    A 2D ray specified by a starting point and an angle
    """

    def __init__(self, z=0, y=0, u=1):
        """
        Args:
            z: axial position of start of ray
            y: distance off axis for start of ray
            u: angle from axial direction
        """
        self.z = z
        self.u = u

    def __str__(self):
        return "z=%.3f y=%.3f u=%.3f" % (self.z, self.y, self.u)

    def __repr__(self):
         return "Ray(z=" + self.z + ", y=" + self.y + ", u=" + self.u + ")"

    
class ParaxialElement:
    """
    A 2D optical element
    """

    def __init__(self, z=0, y=1):
        """
        Args:
            z: axial position
        """
        self.z = z
        self.y = y

    def __str__(self):
        return "z=%.3f" % (self.z)

    def __repr__(self):
         return "Ray(z=" + self.z + ")"
         
    def draw():
        plt.plot([z,z],[-y,y]);


class Mirror(ParaxialElement):
    """
    A mirror element
    """

    def __init__(self, z=0, f=1, y=1):
        """
        Args:
            z: axial position of start of ray
            f: focal length of mirror
        """
        self.z = z
        self.f = f
        self.y = y

    def __str__(self):
        return "z=%.3f y=%.3f f=%.3f" % (self.z, self.y, self.f)

    def __repr__(self):
         return "Mirror(z=" + self.z + ", y=" + self.y + ", f=" + self.f + ")"


class Lens(ParaxialElement):
    """
    A thin lens element
    """

    def __init__(self, z=0, f=1, y=1):
        """
        Args:
            z: axial position of start of ray
            f: focal length of mirror
        """
        self.z = z
        self.f = f
        self.y = y

    def __str__(self):
        return "z=%.3f y=%.3f f=%.3f" % (self.z, self.y, self.f)

    def __repr__(self):
         return "Lens(z=" + self.z + ", y=" + self.y + ", f=" + self.f + ")"

class Plane(ParaxialElement):
    """
    A plane element
    """

    def __init__(self, z=0, y=1):
        """
        Args:
            z: axial position of start of ray
            y: lateral extent of plane
        """
        self.z = z
        self.y = y

    def __str__(self):
        return "z=%.3f y=%.3f" % (self.z, self.y)

    def __repr__(self):
         return "Plane(z=" + self.z + ", y=" + self.y + ")"

class Aperture(Plane):
    """
    A thin lens element
    """

    def __init__(self, z=0, y=1, D=1):
        """
        Args:
            z: axial position of start of ray
            D: diameter of aperture
            y: lateral extent of aperture
        """
        self.z = z
        self.D = D
        self.y = y

    def __str__(self):
        return "z=%.3f D=%.3f y=%.3f" % (self.z, self.D, self.y)

    def __repr__(self):
         return "Lens(z=" + self.z + ", D=" + self.D + ", y=" + self.y + ")"

class Object(Plane):
    """
    An object with height h
    """

    def __init__(self, z=0, h=1):
        """
        Args:
            z: axial position of start of ray
            h: lateral extent of aperture
        """
        self.z = z
        self.y = h

    def __str__(self):
        return "z=%.3f h=%.3f" % (self.z, self.y)

    def __repr__(self):
         return "Object(z=" + self.z + ", D=" + self.D + ", y=" + self.y + ")"

class Image(Plane):
z
h
matrix
draw

class ParaxialSystem:
    """
    A system of optical elements
    """

    def __init__(self, z, elems):
        """
        Args:
            z: array of z-positions of elements
            e: array of optical elements
        """
        self.z = np.empty(0, dtype=float)
        self.e = np.empty(0, dtype=ParaxialElement)
        self.m = np.empty(0, dtype=matrix)
        
    def append(self, z, elem):
        np.append(self.z,z)
        np.append(self.e,elem)
        np.append(self.matrix, elem.matrix)
        
    def move(ray):
        for e in enumerate(self.e):
            e.draw()
            
        for i,e in enumerate(self.z):
            draw_ray_to_next_element
            bend_ray
            
    def draw():
        for e in self.e:
            e.draw()

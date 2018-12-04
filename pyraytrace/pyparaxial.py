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

__all__ = ['Ray',
           'ParaxialElement',
           'Mirror',
           'ThinLens',
           'Plane',
           'Aperture',
           'Object',
           'ParaxialSystem']


class Ray:
    """
    A 2D ray specified by a starting point and an angle
    """

    def __init__(self, z=0, y=0, nu=1):
        """
        Args:
            z: axial position of start of ray
            y: distance off axis for start of ray
            nu: angle from axial direction (radians)
        """
        self.z = z
        self.ynu = np.array([y,nu])

    def __str__(self):
        return "Ray z=%.3f y=%.3f nu=%.3f" % (self.z, self.ynu[0], self.ynu[1])

    def __repr__(self):
         return "Ray(z=" + self.z + ", y=" + self.ynu[0] + ", nu=" + self.ynu[1] + ")"


class ParaxialElement:
    """
    A 2D optical element
    """

    def __init__(self, z=0, h=1):
        """
        Args:
            z: axial position
        """
        self.z = z
        self.h = h
        self.color = 'black'
        self.ABCD = np.array([[1,0], [0,1]])

    def __str__(self):
        return "ParaxialElement z=%.3f, h=%.3f" % (self.z, self.h)

    def __repr__(self):
         return "ParaxialElement(z=" + self.z + ", h=" + self.h + ")"
         
    def draw(self):
        plt.plot([self.z, self.z],[-self.h, self.h],color=self.color)


class Mirror(ParaxialElement):
    """
    A mirror element
    """

    def __init__(self, z=0, h=1, R=np.inf, color='purple'):
        """
        Args:
            z: axial position of start of ray
            R: radius of curvature of mirror
            h: height of mirror (half the diameter)
        """
        self.z = z
        self.R = R
        self.h = h
        self.color = color
        self.ABCD = np.array([[1,0], [2/R, 1]])

    def __str__(self):
        return "Mirror z=%.3f h=%.3f R=%.3f" % (self.z, self.h, self.R)

    def __repr__(self):
        return "Mirror(z=" + self.z + ", h=" + self.h + ", R=" + self.R + ")"

    def draw(self):
        plt.plot([self.z, self.z],[-self.h, self.h],color=self.color)

class ThinLens(ParaxialElement):
    """
    A thin lens element
    """

    def __init__(self, z=0, h=1, f=10, color='blue'):
        """
        Args:
            z: axial position of start of ray
            f: focal length of mirror
            h: height of lens (half the diameter)
        """
        self.z = z
        self.f = f
        self.h = h
        self.color = color
        self.ABCD = np.array([[1,0], [-1/f, 1]])

    def __str__(self):
        return "ThinLens z=%.3f h=%.3f f=%.3f" % (self.z, self.h, self.f)

    def __repr__(self):
         return "ThinLens(z=" + self.z + ", h=" + self.h + ", f=" + self.f + ")"

    def draw(self):
        plt.plot([self.z, self.z],[-self.h, self.h],color=self.color)


class Plane(ParaxialElement):
    """
    A plane element
    """

    def __init__(self, z=0, h=1):
        """
        Args:
            z: axial position of start of ray
            h: lateral extent of plane above optical axis
        """
        self.z = z
        self.h = h
        self.color = 'orange'
        self.ABCD = np.array([[1,0], [0,1]])

    def __str__(self):
        return "Plane z=%.3f h=%.3f" % (self.z, self.h)

    def __repr__(self):
         return "Plane(z=" + self.z + ", h=" + self.h + ")"

    def draw(self):
        plt.plot([self.z, self.z],[-self.h, self.h],color=self.color)


class Aperture(Plane):
    """
    An aperture element
    """

    def __init__(self, z=0, h=1, D=1, width=0.3):
        """
        Args:
            z: axial position of start of ray
            D: diameter of aperture
            h: lateral extent of aperture
            width: fraction of aperture for horizontal bars
        """
        self.z = z
        self.D = D
        self.h = h
        self.width = width
        self.color = 'black'
        self.ABCD = np.array([[1,0], [0,1]])

    def __str__(self):
        return "Aperture z=%.3f h=%.3f D=%.3f width=%.3f" % (self.z, self.h, self.D, self.width)

    def __repr__(self):
         return "Aperture(z=" + self.z + ", h=" + self.h + ", D=" + self.D + ")"

    def draw(self):
        r = self.D/2
        w = self.width * r / 2
        plt.plot([self.z,  self.z  ],[ r, self.h],color=self.color)
        plt.plot([self.z-w,self.z+w],[ r,r],color=self.color)
        plt.plot([self.z,  self.z  ],[-r,-self.h],color=self.color)
        plt.plot([self.z-w,self.z+w],[-r,-r],color=self.color)


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
        self.h = h
        self.color = 'red'
        self.ABCD = np.array([[1,0], [0,1]])

    def __str__(self):
        return "Object z=%.3f h=%.3f" % (self.z, self.h)

    def __repr__(self):
         return "Object(z=" + self.z + ", h=" + self.h + ")"

    def draw(self):
        plt.plot([self.z, self.z],[0,self.h],color=self.color)
        

class ParaxialSystem:
    """
    A system of optical elements
    """

    def __init__(self):
        self.elements = []
        
    def __str__(self):
        str = ''
        for e in self.elements:
            str = str + e.__str__() + "\n"
        return str

    def __repr__(self):
        repr = ''
        for e in self.elements:
            repr = repr + e.__repr__() + "\n"
        return repr

    def append(self, elem):
        self.elements.append(elem)
        
    def move(self, ray):
        self.draw()
            
        for i,e in enumerate(self.elements):
            print(i,e)
            
    def draw(self):
        plt.plot([self.elements[0].z,self.elements[-1].z],[0,0],'k:')
        for e in self.elements:
            e.draw()

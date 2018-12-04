# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
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

    def __init__(self, y=0, nu=1):
        """
        Args:
            y: distance off axis for start of ray
            nu: angle from axial direction (radians)
        """
        self.ynu = np.array([y, nu], dtype=float)

    def __str__(self):
        return "Ray y=%.3f nu=%.3f" % (self.ynu[0], self.ynu[1])

    def __repr__(self):
        return "Ray(y=%f, nu=%f)" % (self.ynu[0], self.ynu[1])

    def aim(self, z0=0, y0=0, z1=1, y1=1):
        """ aim ray from (z0,y0) to (z1,q1) """
        self.ynu[0] = y0
        self.ynu[1] = (y1 - y0) / (z1 - z0)

    def aim2(self, p, q):
        """ aim ray from point p to point q """
        self.ynu[0] = p[1]
        self.ynu[1] = (q[1] - p[1]) / (q[0] - p[0])


class ParaxialElement:
    """
    A 2D optical element
    """

    def __init__(self, z=0, h=1, color='black'):
        """
        Args:
            z: axial position
        """
        self.z = z
        self.h = h
        self.ABCD = np.array([[1, 0], [0, 1]], dtype=float)
        self.color = color

    def __str__(self):
        return "ParaxialElement z=%.3f, h=%.3f" % (self.z, self.h)

    def __repr__(self):
        return "ParaxialElement(z=%f, h=%f)" % (self.z, self.h)

    def draw(self):
        """ Draw the element """
        plt.plot([self.z, self.z], [-self.h, self.h], color=self.color)


class Mirror(ParaxialElement):
    """
    A mirror element
    """

    def __init__(self, z=0, h=1, R=np.inf, color='purple'):
        """
        Args:
            z: axial position of mirror
            R: radius of curvature of mirror
            h: height of mirror (half the diameter)
        """
        super().__init__(z=z, h=h, color=color)
        self.R = R
        self.ABCD = np.array([[1, 0], [2/R, 1]], dtype=float)

    def __str__(self):
        return "Mirror z=%.3f h=%.3f R=%.3f" % (self.z, self.h, self.R)

    def __repr__(self):
        return "Mirror(z=%f, h=%f, R=%f)" % (self.z, self.h, self.R)

    def draw(self):
        """ Draw the mirror """
        plt.plot([self.z, self.z], [-self.h, self.h], color=self.color)


class ThinLens(ParaxialElement):
    """
    A thin lens element
    """

    def __init__(self, z=0, h=1, f=10, width=0.2, color='blue'):
        """
        Args:
            z: axial position of thin lens
            f: focal length of mirror
            h: height of lens (half the diameter)
            width: width of arrowhead
        """
        super().__init__(z=z, h=h, color=color)
        self.f = f
        self.width = width
        self.ABCD = np.array([[1, 0], [-1/f, 1]], dtype=float)

    def __str__(self):
        return "ThinLens z=%.3f h=%.3f f=%.3f" % (self.z, self.h, self.f)

    def __repr__(self):
        return "ThinLens(z=%f, h=%f, f=%f)" % (self.z, self.h, self.f)

    def draw(self):
        """ Draw the thin lens """
        w = self.h * self.width
        plt.plot([self.z, self.z], [-self.h, self.h], color=self.color)
        if self.f > 0:
            plt.plot([self.z-w, self.z, self.z+w],
                     [self.h-w, self.h, self.h-w], color=self.color)
            plt.plot([self.z-w, self.z, self.z+w],
                     [-self.h+w, -self.h, -self.h+w], color=self.color)
        else:
            plt.plot([self.z-w, self.z, self.z+w],
                     [self.h+w, self.h, self.h+w], color=self.color)
            plt.plot([self.z-w, self.z, self.z+w],
                     [-self.h-w, -self.h, -self.h-w], color=self.color)


class Plane(ParaxialElement):
    """
    A plane element
    """

    def __init__(self, z=0, h=1, color='black'):
        """
        Args:
            z: axial position of plane
            h: lateral extent of plane above optical axis
        """
        super().__init__(z=z, h=h, color=color)

    def __str__(self):
        return "Plane z=%.3f h=%.3f" % (self.z, self.h)

    def __repr__(self):
        return "Plane(z=%f, h=%f)" % (self.z, self.h)

    def draw(self):
        """ Draw the plane """
        plt.plot([self.z, self.z], [-self.h, self.h], ':', color=self.color)


class Aperture(Plane):
    """
    An aperture element
    """

    def __init__(self, z=0, h=1, D=1, width=0.3, color='black'):
        """
        Args:
            z: axial position of aperture
            D: diameter of aperture
            h: overall height of aperture (needed for drawing)
            width: fraction of aperture for horizontal bars
        """
        super().__init__(z=z, h=h, color=color)
        self.D = D
        self.width = width

    def __str__(self):
        return "Aperture z=%.3f h=%.3f D=%.3f width=%.3f" % (self.z, self.h, self.D, self.width)

    def __repr__(self):
        return "Aperture(z=%f, h=%f, D=%f)" % (self.z, self.h, self.D)

    def draw(self):
        r = self.D/2
        w = self.width * self.h / 2
        plt.plot([self.z, self.z], [r, self.h], color=self.color)
        plt.plot([self.z-w, self.z+w], [r, r], color=self.color)
        plt.plot([self.z, self.z], [-r, -self.h], color=self.color)
        plt.plot([self.z-w, self.z+w], [-r, -r], color=self.color)


class Object(Plane):
    """
    An object with height h
    """

    def __init__(self, z=0, h=1, color='red'):
        """
        Args:
            z: axial position of start of ray
            h: lateral extent of aperture
        """
        super().__init__(z=z, h=h, color=color)

    def __str__(self):
        return "Object z=%.3f h=%.3f" % (self.z, self.h)

    def __repr__(self):
        return "Object(z=%f, h=%f)" % (self.z, self.h)

    def draw(self):
        #        plt.plot([self.z, self.z], [0, self.h], color=self.color)
        plt.annotate("", xy=(self.z, 0), xytext=(self.z, self.h),
                     arrowprops=dict(arrowstyle="<-"))


class ParaxialSystem:
    """
    A system of optical elements
    """

    def __init__(self):
        self.elements = []

    def __str__(self):
        s = ''
        for e in self.elements:
            s += e.__str__() + "\n"
        return s

    def __repr__(self):
        r = ''
        for e in self.elements:
            r += e.__repr__() + "\n"
        return r

    def append(self, elem):
        """ add another optical element to the system """
        self.elements.append(elem)

    def minmax(self):
        """ return min and max z """
        zmin = 1e6
        zmax = -1e6
        for e in self.elements:
            if e.z < zmin:
                zmin = e.z
            if zmax < e.z:
                zmax = e.z
        return zmin, zmax

    def draw_ray_through_system(self, ray, color='black'):
        """ Draw the ray as it moves through the system """

        z_last = None
        ynu = np.copy(ray.ynu)

        for e in self.elements:
            if z_last != None:
                prop = np.array([[1, e.z-z_last], [0, 1]])

                ynu_next = np.matmul(prop, ynu)
                plt.plot([z_last, e.z], [ynu[0], ynu_next[0]], color=color)

                ynu_bent = np.matmul(e.ABCD, ynu_next)
                ynu = np.copy(ynu_bent)
            z_last = e.z

    def draw(self):
        """ Draw all the elements """
        zmin, zmax = self.minmax()

        plt.plot([zmin, zmax], [0, 0], 'k')
        for e in self.elements:
            e.draw()

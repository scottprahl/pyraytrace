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

class Plane:
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
        a = "xyz=[%.3f,%.3f,%.3f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "uvw=[%.3f,%.3f,%.3f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        length = np.dot(self.uvw, self.uvw)
        return a + ", " + b + " norm=%.4f" % length + " D=%f" % self.D

    def __repr__(self):
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "[%f,%f,%f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
        return "Plane(" + a + ", " + b + ")"

    def draw_zy(self, ymin=0, ymax=1, zmin=0, zmax=1):
        """
        Draw representation in the zy-plane (x==0) that lies in
        the rectangle bounded by ymin,ymax and zmin,zmax

        Thus  v*y + w*z + D = 0

        """
        if self.uvw[2] != 0:
            ymn = ymin
            ymx = ymax
            zmn = -(self.D + self.uvw[1] * ymin) / self.uvw[2]
            zmx = -(self.D + self.uvw[1] * ymax) / self.uvw[2]         
            zmn = max(zmin,zmn)
            zmx = min(zmax,zmx)
            print("   zy=(%.2f,%.2f), zy=(%.2f,%.2f)"%(zmn,ymn,zmx,ymx))
            plt.plot([zmn, zmx], [ymn, ymx], 'k')
            return
            
        if self.uvw[1] != 0:
            ymn = -(self.D + self.uvw[2] * zmin) / self.uvw[1]
            ymx = -(self.D + self.uvw[2] * zmax) / self.uvw[1]
            ymn = max(ymn,ymin)
            ymx = min(ymx,ymax)
            zmn = zmin
            zmx = zmax
            print("   zy=(%.2f,%.2f), zy=(%.2f,%.2f)"%(zmn,ymn,zmx,ymx))
            plt.plot([zmn, zmx], [ymn, ymx], 'k')


    def distance(self, ray):
        """
        distance from start of ray to plane
        """
        cos_angle = np.dot(ray.uvw, self.uvw)
        if abs(cos_angle) < 1e-8:
            return np.inf

        return -(np.dot(ray.xyz, self.uvw) + self.D) / cos_angle

    def is_in_plane(self, point):
        """
        return True/False if point is in the plane
        """
        dist = abs(np.dot(point, self.uvw) + self.D)
        return dist < 1e-6


class Ray:
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
        a = "xyz=[%.3f,%.3f,%.3f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = "uvw=[%.3f,%.3f,%.3f]" % (self.uvw[0], self.uvw[1], self.uvw[2])
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

def refract(uvw, normal, ni, nt):
        """
        Spencer and Murty, equation 36
        """
        cosine = np.dot(normal, uvw)
        if cosine < 0:
            cosine *= -1
            normal *= -1
        
        refractive_index = nt/ni
        a = refractive_index * cosine
        b = refractive_index ** 2 - 1
        disc = a ** 2 - b
        if disc < 0:  # reflected
            out = uvw - 2 * cosine * normal
        else:
            g = -a + np.sqrt(disc)
            out = refractive_index * uvw + g * normal

        return out

class Sphere:
    """
    A class to help to ray-trace spherical objects

    A sphere is defined by the equation:

        (x-x0)**2 + (y-y0)**2 + (z-z0)**2 = R**2

    where (x0,y0,z0) is the center of the sphere and R is the radius
    """

    def __init__(self, xyz=(0, 0, 0), R=1.0, n= 1.0):
        """
        Args:
            xyz: array describing the center of the sphere in cartesian coordinates
            R:   radius of the sphere
            n:   index of refraction of the sphere
        """

        self.xyz = np.array(xyz)
        self.radius = R
        self.n = n

    def __str__(self):
        a = "center=[%.3f,%.3f,%.3f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        b = ", radius = %f" % self.radius
        return a + b

    def __repr__(self):
        a = "[%f,%f,%f]" % (self.xyz[0], self.xyz[1], self.xyz[2])
        return "Sphere(" + a + ", %f" % self.radius + ")"

    def draw_zy(self, ymax=np.inf, side='both'):
        """
        Draw representation in the zy-plane
        """
        RR = np.sqrt(self.radius**2 - self.xyz[0]**2)
        yy = min(ymax, RR)
        y = np.linspace(-yy, yy, 50)
        r = RR**2 - (y - self.xyz[1])**2
        np.place(r, r < 0, 0)
        z = np.sqrt(r)
        if side == 'both' or side == 'right':
            plt.plot(z + self.xyz[2], y, 'k')
        if side == 'both' or side == 'left':
            plt.plot(-z + self.xyz[2], y, 'k')

    def unit_normal_at(self, point):
        """
        Return outward normal to point on sphere
        """
        diff = point - self.xyz
        mag = np.sqrt(np.dot(diff, diff))
        return diff / mag

    def distance(self, ray):
        """
        Return the nearest positive distance of a ray to the sphere
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
        """
        Spencer and Murty, equation 36
        """
        normal = self.unit_normal_at(ray.xyz)
        if outside:
            return refract(ray.uvw, normal, 1, n)
        else:
            return refract(ray.uvw, normal, n, 1)


class Prism:
    """
    A class to help to ray-trace through prisms

    A prism is defined by three planes
    """

    def __init__(self, A, B, C):
        """
        Args:
            A: plane object for first side
            B: plane object for second side
            C: plane object for third side
        """
        self.A = A
        self.B = B
        self.C = C

    def __str__(self):
        return self.A.__str__() + '\n' + self.B.__str__() + '\n' + self.C.__str__()

    def __repr__(self):
        return self.A.__repl__() + self.B.__repl__() + self.C.__repl__()

    def draw_zy(self, ymin=0, ymax=1, zmin=0, zmax=1):
        """
        Draw representation in the zy-plane

        each plane satisfies u*x + v*y + w*z + D = 0.  In the zy-plane x==0
        therefore

        v*y + w*z + D = 0
        the corners (y,z) can be found by solving

        v1*y + w1*z + D1 = 0
        v2*y + w2*z + D2 = 0

        y = -(D1*w2-D2*w1)/(w2*v1-w1*v2)
        z =  (D1*v2-D2*v1)/(w2*v1-w1*v2)
        """

        denom = self.B.uvw[2]*self.A.uvw[1]-self.A.uvw[2]*self.B.uvw[1]
        y1 = -(self.A.D*self.B.uvw[2]-self.B.D*self.A.uvw[2])/denom
        z1 = -(self.A.D*self.B.uvw[1]-self.B.D*self.A.uvw[1])/denom

        self.A.draw_zy(ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
        self.B.draw_zy(ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
        self.C.draw_zy(ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)

    def unit_normal_at(self, point):
        """
        Return outward normal to point on prism
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
        """
        Return the nearest positive distance of a ray to a prism face
        """
        d1 = self.A.distance(ray)
        d2 = self.B.distance(ray)
        d3 = self.C.distance(ray)
        dd = np.array([d1, d2, d3])
        np.place(dd, dd <= 1e-6, 999)
        # print("side 1 d=%.3f\n"%dd[0])
        # print("side 2 d=%.3f\n"%dd[1])
        # print("side 3 d=%.3f\n"%dd[2])
        return min(dd)

    def refract(self, ray, refractive_index):
        """
        Spencer and Murty, equation 36
        """
        normal = self.unit_normal_at(ray.xyz)
        if outside:
            return refract(ray.uvw, normal, 1, n)
        else:
            return refract(ray.uvw, normal, n, 1)


class Lens:
    """
    A class to help to ray-trace through a lens

    A lens is defined by two surfaces
    """

    def __init__(self, surface1, surface2, refractive_index, thickness):
        """
        Args:
            surface1: first surface
            surface2: second surface
            refractive_index:        index of refraction
            d:        thickness of lens
        """
        self.surface1 = surface1
        self.surface2 = surface2
        self.refractive_index = refractive_index
        self.thickness = thickness

    def __str__(self):
        a = str(self.surface1)
        b = str(self.surface2)
        c = "refractive index = %f" % self.refractive_index
        d = "thickness = %f" % self.thickness
        return a + "\n" + b + "\n" + c + "\n" + d

    def __repr__(self):
        a = repr(self.surface1)
        b = repr(self.surface2)
        c = ", %f, %f)" % (self.refractive_index, self.thickness)
        return "Lens(" + a + "," + b + c

    def distance(self, ray, which_surface):
        """
        Distance to surface
        """
        if which_surface == 1:
            return self.surface1.distance(ray)

        return self.surface2.distance(ray)

    def refract(self, ray, which_surface):
        """
        Bend light at surface
        """
        if which_surface == 1:
            return self.surface1.refract(ray, 1/self.refractive_index)

        return self.surface2.refract(ray, self.refractive_index)

    def draw_zy(self):
        """
        Draw representation in the zy-plane
        """
        self.surface1.draw_zy(side='left')
        self.surface2.draw_zy(side='right')


class ThinLens:
    """
    A class for a thin lens
    """
    def __init__(self, focal_length, vertex, diameter=10):
        """
        Args:
            vertex: first surface
            f: focal length
            diameter: diameter of lens
        """
        self.f = focal_length
        self.vertex = vertex
        self.diameter = diameter

    def __str__(self):
        a = "focal length = %f" % self.f
        b = "vertex = %f" % self.vertex
        c = "diameter = %f" % self.diameter
        return a + "\n" + b + "\n" + c

    def __repr__(self):
        return "ThinLens(%f, %f, diameter=%f)" % (self.f, self.vertex, self.diameter)
        
    def distance(self, ray):
        """
        Distance to the lens
        """
        z0 = ray.xyz[2]
        z1 = self.vertex
        w = ray.uvw[2]
        
        if w==0:
            return np.inf
        else :
            return (z1-z0)/w
        
    def refract(self, ray):
        """
        Bend light at surface
        """
        w = ray.uvw[2]
        yf = self.f/w
        y0 = ray.xyz[1]
        return self.f/np.sqrt((yf-y0)**2+self.f**2)

    def draw_zy(self):
        """
        Draw representation in the zy-plane
        """
        plt.plot([self.vertex,self.vertex], [-self.diameter/2,self.diameter/2], ':r')
  
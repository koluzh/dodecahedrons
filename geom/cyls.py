from geom.core import *
import gmsh
import time
import numpy as np


class Cylinder:
    def __init__(self, center: Dot, r: float, h: float, n: int, angle: Dot = None):
        self.center = center
        self.r = r
        self.h = h
        self.n = n
        if angle is None:
            angle = Dot(0, 0, 0)
        self.angle = angle
        rot = get_rotation_matrix(angle)
        v = np.matrix([1, 0, 0]) @ rot
        v = np.ravel(v)
        self.v = Dot(coords=v)

    def build_mesh(self):
        x, y, z = self.center.as_list()
        dx, dy, dz = Dot(coords=self.v.coords * self.h).as_list()
        gmsh.model.occ.add_cylinder(x, y, z, dx, dy, dz, self.r, self.n)


def cyl_intersection(c1: Cylinder, c2: Cylinder) -> bool:
    pass
from main import Dodecahedron
from main import Dot
import numpy as np
from numba import njit


class Pentagon:
    def __init__(self, dots):
        self.dots = dots
        self.vectors = []
        self.build_vectors()
        self.plane = None
        self.plane_D = None
        self.build_plane()
        self.check_plane()

    def build_vectors(self):
        for i in range(0, 4):
            temp_v = self.dots[i + 1].coords - self.dots[i].coords
            self.vectors.append(temp_v)
        temp_v = self.dots[0].coords - self.dots[4].coords
        self.vectors.append(temp_v)

    def build_plane(self):
        v1 = self.vectors[0]
        v2 = self.vectors[1]
        plane = np.cross(np.squeeze(v1), np.squeeze(v2))
        d = np.dot(plane, self.dots[0].coords)
        self.plane = plane
        self.plane_D = d

    def check_plane(self):
        for d in self.dots:
            temp_val = np.dot(d.coords, self.plane)
            if temp_val != d:
                Exception("dots are not in the same plane")

    def check_dot_on_plane(self, d: Dot):
        temp_val = np.dot(d.coords, self.plane)
        if temp_val == self.D:
            return True
        else:
            return False

    def check_line_intersection(self, d0: Dot, d1: Dot, epsilon=1e-6):
        p0 = d0.coords
        p1 = d1.coords
        u = p1 - p0
        dot = np.dot(self.plane, u)

        if abs(dot) > epsilon:
            p_co = self.plane * (-self.plane_D / np.linalg.norm(self.plane))
            w = p0 - p_co
            fac = -np.dot(self.plane, w) / dot
            u = u * fac
            return p0 + u

        return None

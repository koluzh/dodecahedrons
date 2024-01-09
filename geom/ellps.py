import numpy as np
import numpy.linalg

from geom.core import *
import time


class Ellipsoid:
    # kx^2 + ly^2 + mz^2 - r^2 = 0
    # n is number of ellipsoid needed for gmsh
    def __init__(self, center: Dot, k: float, l: float, m: float, r: float, angle: Dot, num: int):
        self.center = center
        # radius
        self.r = r
        # describes ellipsoid
        self.k = k
        self.l = l
        self.m = m
        # angles
        self.angle = angle
        self.alpha = angle.x
        self.beta = angle.y
        self.gamma = angle.z
        # number of ellipsoid
        self.num = num

        self.a = r / np.sqrt(k)
        self.b = r / np.sqrt(l)
        self.c = r / np.sqrt(m)
        self.volume = self.a * self.b * self.c * np.pi * 4 / 3
        klm = [self.k, self.l, self.m]
        self.info = [self.center.coords, klm, self.angle.coords, self.r]
        self.matrix = self.get_matrix()

    def create_mesh(self):
        gmsh.model.occ.add_sphere(self.center.x, self.center.y, self.center.z, self.r, tag=self.num)
        gmsh.model.occ.dilate([(3, self.num)], self.center.x, self.center.y, self.center.z, 1 / np.sqrt(self.k),
                              1 / np.sqrt(self.l), 1 / np.sqrt(self.m))
        gmsh.model.occ.mesh.set_size(gmsh.model.occ.get_entities(0), 0.5)
        gmsh.model.occ.rotate([(3, self.num)], self.center.x, self.center.y, self.center.z,
                              1, 0, 0, self.alpha)
        gmsh.model.occ.rotate([(3, self.num)], self.center.x, self.center.y, self.center.z,
                              0, 1, 0, self.beta)
        gmsh.model.occ.rotate([(3, self.num)], self.center.x, self.center.y, self.center.z,
                              0, 0, 1, self.gamma)

        # yay

    def in_box(self, box: Box, epsilon: float = None):
        if epsilon is None:
            epsilon = 0
        d1 = box.s
        d2 = box.d
        for i in range(3):
            if self.center.coords[i] <= d1.coords[i] - epsilon or self.center.coords[i] >= d2.coords[i] + epsilon:
                return False

        for p in box.planes:
            temp_dot = p.get_projection(self.center)
            # gmsh.model.occ.add_point(temp_dot.x, temp_dot.y, temp_dot.z)
            # print(temp_dot.coords)
            t = d_in_ell(temp_dot, self)
            if t > 0:
                t = np.sqrt(t)
            if t <= epsilon:
                return False
        return True

    def print_info(self):
        print(self.num, 'th ell center, klm, angle, r')
        for i in self.info:
            print(i)

    def get_matrix(self) -> np.matrix:
        A = self.c ** 2 * self.b ** 2
        B = self.a ** 2 * self.c ** 2
        C = self.a ** 2 * self.b ** 2
        K = -1 * self.a ** 2 * self.b ** 2 * self.c ** 2
        init_mat = np.matrix(f"{A} 0 0 0; 0 {B} 0 0; 0 0 {C} 0; 0 0 0 {K}")
        rotated_matrix = init_mat @ get_rotation_matrix_4(self.angle)
        translation_matrix = np.matrix(f"1 0 0 {-self.center.x}; 0 1 0 {-self.center.y}; 0 0 1 {-self.center.z}; 0 0 0 1")
        # print(init_mat, '\n')
        # print(f"ell {self.num}")
        # print(get_rotation_matrix_4(self.angle), '\n')
        # print(translation_matrix, '\n')
        # print(translation_matrix.T, '\n')
        TRS = translation_matrix @ rotated_matrix @ translation_matrix.T
        return TRS


def get_rotation_matrix_4(d: Dot) -> np.matrix:
    r_mat = get_rotation_matrix(d)
    r = r_mat
    temp1 = np.array(r[0, :])
    temp1 = temp1.ravel()
    temp1 = np.hstack([temp1, [0]])
    temp2 = np.array(r[1, :])
    temp2 = temp2.ravel()
    temp2 = np.hstack([temp2, [0]])
    temp3 = np.array(r[2, :])
    temp3 = temp3.ravel()
    temp3 = np.hstack([temp3, [0]])
    r = np.matrix([temp1, temp2, temp3, [0, 0, 0, 1]])
    # print(r)
    r_mat = r
    return r_mat

def rot_dot(d: Dot, angle: Dot, reverse: bool = None) -> Dot:
    d2_t = np.array([d.x, d.y, d.z, 1])
    r_mat = get_rotation_matrix_4(angle)
    d2_t = d2_t @ r_mat
    d2_t = np.array([d2_t[0, 0], d2_t[0, 1], d2_t[0, 2]])  # this is so shit
    d2_t = Dot(coords=d2_t)
    return d2_t


def d_to_ecs(d: Dot, e: Ellipsoid, reverse: bool = None):
    # Dot coordinates translated to ellipsoid coordinate system
    c = e.center
    d_c_t = np.array(d.coords) - np.array(c.coords)
    d_t = Dot(coords=d_c_t)
    d_t_r = rot_dot(d_t, e.angle, reverse)
    return d_t_r


def d_in_ell(d: Dot, e: Ellipsoid):
    d = d_to_ecs(d, e)
    x, y, z = d.x, d.y, d.z
    t = e.k * (x ** 2) + e.l * (y ** 2) + e.m * (z ** 2) - e.r ** 2
    return t

def ell_2_ell(e1: Ellipsoid, e2: Ellipsoid) -> bool:
    pass
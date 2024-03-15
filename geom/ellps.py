import gmsh
import numpy as np
import numpy.linalg

from geom.core import *
import time

from scipy.linalg import eigh
from scipy.optimize import minimize_scalar


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

    def build_mesh(self):
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

    def ellipsoid_intersection_test_helper(self,
                                    Sigma_A: np.matrix,
                                    Sigma_B: np.ndarray,
                                    mu_A: np.ndarray,
                                    mu_B: np.ndarray):
        lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
        v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
        return lambdas, Phi, v_squared

    def ellipsoid_K_function(self, ss, lambdas, v_squared, tau):
        ss = np.array(ss).reshape((-1, 1))
        lambdas = np.array(lambdas).reshape((1, -1))
        v_squared = np.array(v_squared).reshape((1, -1))
        return 1. - (1. / tau ** 2) * np.sum(v_squared * ((ss * (1. - ss)) / (1. + ss * (lambdas - 1.))), axis=1)

    def ellipsoid_intersection_test(self,
                                    Sigma_A: np.matrix,
                                    Sigma_B: np.matrix,
                                    mu_A: np.ndarray,
                                    mu_B: np.ndarray,
                                    tau: float) -> bool:
        eigvals, eigvecs = np.linalg.eig(Sigma_B)
        l = np.diag(eigvals)
        # print(l)
        lambdas, Phi, v_squared = self.ellipsoid_intersection_test_helper(Sigma_A, l, mu_A, mu_B)
        res = minimize_scalar(self.ellipsoid_K_function,
                              bracket=[0.0, 0.5, 1.0],
                              args=(lambdas, v_squared, tau))
        return res.fun[0] >= 0 - 0.000_001

    def intersects(self, ell: 'Ellipsoid') -> bool:
        return self.ellipsoid_intersection_test(
            np.linalg.inv(self.matrix),
            np.linalg.inv(ell.matrix),
            self.center.coords,
            ell.center.coords, tau=1)

    def print_info(self):
        print(self.num, 'th ell center, klm, angle, r')
        for i in self.info:
            print(i)

    def get_matrix(self) -> np.matrix:
        A = 1 / self.a ** 2
        B = 1 / self.b ** 2
        C = 1 / self.c ** 2
        init_mat = np.diag([A, B, C])
        TRS = get_rotation_matrix(self.angle) @ init_mat @ get_rotation_matrix(self.angle).T
        return TRS

    def in_box(self, box: Box) -> bool:
        # translate dots of plane to ecs, translate ecs so that ellipsoid becomes sphere,
        # check if sphere collides with plane
        planes = box.planes
        planes_in_ecs = [plane_to_ecs(p, self) for p in planes]

        squeeze_mat = np.diag([1 / self.a, 1 / self.b, 1 / self.c])

        planes_squeezed = [plane_squeeze(p, squeeze_mat) for p in planes_in_ecs]
        doesnt_intersect = [not self.intersects_plane(p) for p in planes_squeezed]
        return all(doesnt_intersect)

    def intersects_plane(self, plane_translated: Plane):
        denom = plane_translated.x ** 2 + plane_translated.y ** 2 + plane_translated.z ** 2
        denom = np.sqrt(denom)
        dist = abs(plane_translated.d / denom)
        # print(dist)
        if dist <= self.r:
            return True
        return False


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
    r_mat = r
    return r_mat


def rot_dot(d: Dot, angle: Dot, reverse: bool = None) -> Dot:
    d2_t = np.array([d.x, d.y, d.z, 1])
    r_mat = get_rotation_matrix_4(angle)
    d2_t = d2_t @ np.linalg.inv(r_mat)
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


def d_in_ell(d: Dot, e: Ellipsoid) -> bool:
    d = d_to_ecs(d, e)
    x, y, z = d.x, d.y, d.z
    t = e.k * (x ** 2) + e.l * (y ** 2) + e.m * (z ** 2) - e.r ** 2
    if t <= 1:
        result = True
    else:
        result = False
    return result


def plane_to_ecs(p: Plane, e: Ellipsoid) -> Plane:
    dots = p.dots_list
    dots_traversed = [Dot(coords=d.coords - e.center.coords) for d in dots]
    dots_in_ecs = [d_to_ecs(d, e) for d in dots_traversed]
    return Plane(dots_in_ecs)


def plane_squeeze(p: Plane, squeeze_mat: np.ndarray) -> Plane:
    dots = p.dots_list
    dots_squeezed = [Dot(coords=d.coords @ squeeze_mat) for d in dots]
    return Plane(dots_squeezed)

if __name__ == "__main__":
    import sys
    filename = 'GFG.msh'
    gmsh.initialize()
    e = Ellipsoid(Dot(2.47, 4.01, 3.34), k=1.72, l=1.42, m=1.85, r=1.99, num=1, angle=Dot(-0.89, -0.8, 1.46))
    e.build_mesh()
    a = Dot(0, 0, 10)
    a.build_mesh()
    b = Dot(0, 10, 0)
    b.build_mesh()
    c = Dot(0, 10, 10)
    c.build_mesh()
    p = Plane(d=[a, b, c])
    print(e.intersects_plane(p))
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate()
    gmsh.write(filename)

    if 'close' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

import numpy as np
import geom
from geom import Dot
import gmsh
import sys


class Line:
    def __init__(self, p: Dot, d: Dot):
        self.p = p
        # temp_d = d.coords
        # temp_d = temp_d / np.linalg.norm(temp_d)
        # temp_d = Dot(coords=temp_d)
        self.d = d

    def get_dot(self, t: float):
        temp_dot = self.p.coords + self.d.coords * t
        temp_dot = Dot(coords=temp_dot)
        return temp_dot


class Ellipsoid:
    # kx^2 + ly^2 + mz^2 - r^2 = 0
    # n is number of ellipsoid needed for gmsh
    def __init__(self, center: Dot, k: float, l: float, m: float, r: float, angle: Dot, n: int):
        self.center = center
        # radius
        self.r = r
        # describes ellipsoid
        self.k = k
        self.l = l
        self.m = m
        # angles
        self.alpha = angle.x
        self.beta = angle.y
        self.gamma = angle.z
        # number of ellipsoid
        self.n = n
        self.rot_mat = self.get_rotation_matrix()

    def create_mesh(self):
        gmsh.model.occ.add_sphere(self.center.x, self.center.y, self.center.z, self.r, tag=self.n)
        gmsh.model.occ.dilate([(3, self.n)], self.center.x, self.center.y, self.center.z, 1/np.sqrt(self.k),
                              1/np.sqrt(self.l), 1/np.sqrt(self.m))

    def get_rotation_matrix(self):
        r_x = np.matrix([[1, 0, 0],
                         [0, np.cos(self.alpha), -np.sin(self.alpha)],
                         [0, np.sin(self.alpha), np.cos(self.alpha)]])

        r_y = np.matrix([[np.cos(self.beta), 0, np.sin(self.beta)],
                         [0, 1, 0],
                         [-np.sin(self.beta), 0, np.cos(self.beta)]])

        r_z = np.matrix([[np.cos(self.gamma), -np.sin(self.gamma), 0],
                         [np.sin(self.gamma), np.cos(self.gamma), 0],
                         [0, 0, 1]])

        r = r_x @ r_y
        r = r @ r_z
        # r is rotation matrix
        # bs to extend 3x3 matrix to 4x4 matrix
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
        return r
        # yay


def line_ellipsoid_intersection(e: Ellipsoid, l: Line):
    num_sol = 0
    # copied from geometric tools for computer graphics p. 504
    a = e.k * (l.d.x ** 2) + e.l * (l.d.y ** 2) + e.m * (l.d.z ** 2)
    b = 2 * e.k * (l.p.x) * l.d.x + 2 * e.l * (l.p.y) * l.d.y\
        + 2 * e.m * (l.p.z) * l.d.z
    c = e.k * (l.p.x) ** 2 + e.l * (l.p.y) ** 2 + e.m * (l.p.z) ** 2 - e.r ** 2
    print('r', e.r)
    print('abc', a, b, c)
    if a == 0:
        t = -c/b
        return l.get_dot(t)
    dscrm = b ** 2 - a * c * 4
    print('d=',dscrm)

    t = list()

    if dscrm > 0:
        t.append((-b + np.sqrt(dscrm)) / 2 / a)
        t.append((-b - np.sqrt(dscrm)) / 2 / a)
    elif dscrm == 0:
        t.append((-b / 2 / a))

    intersections = list()

    for i in t:
        print('t = ', i)
        intersections.append(l.get_dot(i))
    return intersections


def create_line(d1: Dot, d2: Dot):
    v = d2.coords - d1.coords
    return Line(d1, Dot(coords=v))


gmsh.initialize()

filename = "GFG.msh"

# origin point
gmsh.model.occ.add_point(0, 0, 0, 0.001)
gmsh.model.occ.add_point(0, 2, 0)
gmsh.model.occ.add_point(0, (2 + np.sqrt(2)), 0)

x = 2 ** (1/3)
kek = Ellipsoid(Dot(0, 0, 0), 1, 1, 1, 2, Dot(0, 0, 0), 1)
kek.create_mesh()

lol = Ellipsoid(Dot(-1, (2 + np.sqrt(2)), 0), 2, 1, 1, 2, Dot(0, 0, 0), 2)
lol.create_mesh()

lin = create_line(lol.center, kek.center)
print(lin.p.coords)
print(lin.d.coords)
k = line_ellipsoid_intersection(kek, lin)

for x in k:
    print(x.coords)


gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write(filename)

if 'close' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

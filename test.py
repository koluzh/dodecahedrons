import numpy as np
from geom import Dot
import gmsh
import random
import sys
from scipy.spatial.transform import Rotation as rot


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
        self.angle = angle
        self.alpha = angle.x
        self.beta = angle.y
        self.gamma = angle.z
        # number of ellipsoid
        self.n = n

    def create_mesh(self):
        gmsh.model.occ.add_sphere(self.center.x, self.center.y, self.center.z, self.r, tag=self.n)
        gmsh.model.occ.dilate([(3, self.n)], self.center.x, self.center.y, self.center.z, 1/np.sqrt(self.k),
                              1/np.sqrt(self.l), 1/np.sqrt(self.m))
        gmsh.model.occ.mesh.set_size(gmsh.model.occ.get_entities(0), 0.1)
        gmsh.model.occ.rotate([(3, self.n)], self.center.x, self.center.y, self.center.z,
                            1, 0, 0, self.alpha)
        gmsh.model.occ.rotate([(3, self.n)], self.center.x, self.center.y, self.center.z,
                            0, 1, 0, self.beta)
        gmsh.model.occ.rotate([(3, self.n)], self.center.x, self.center.y, self.center.z,
                            0, 0, 1, self.gamma)

    def in_box(self, box: Box):
        pass




def get_rotation_matrix(angle: Dot):
    r = rot.from_euler('xyz', angle.coords)
    return r.as_matrix()

def line_ellipsoid_intersection(e: Ellipsoid, l: Line):
    num_sol = 0
    # copied from geometric tools for computer graphics p. 504
    a = e.k * (l.d.x ** 2) + e.l * (l.d.y ** 2) + e.m * (l.d.z ** 2)
    b = 2 * e.k * (l.p.x) * l.d.x + 2 * e.l * (l.p.y) * l.d.y\
        + 2 * e.m * (l.p.z) * l.d.z
    c = e.k * (l.p.x) ** 2 + e.l * (l.p.y) ** 2 + e.m * (l.p.z) ** 2 - e.r ** 2
    #print('r', e.r)
    #print('abc', a, b, c)
    if a == 0:
        t = -c/b
        return l.get_dot(t)
    dscrm = b ** 2 - a * c * 4
    # print('d=',dscrm)

    t = list()

    if dscrm > 0:
        t.append((-b + np.sqrt(dscrm)) / 2 / a)
        t.append((-b - np.sqrt(dscrm)) / 2 / a)
    elif dscrm == 0:
        t.append((-b / 2 / a))

    intersections = list()

    for i in t:
        # print('t = ', i)
        intersections.append(l.get_dot(i))

    return intersections


def create_line(d1: Dot, d2: Dot):
    v = d2.coords - d1.coords
    return Line(d1, Dot(coords=v))


def rot_dot(d: Dot, angle: Dot, reverse: bool = None):
    d2_t = np.array([d.x, d.y, d.z, 1])
    r_mat = get_rotation_matrix(angle)
    if reverse is True:
        r_mat = np.linalg.inv(r_mat)
    # bs to extend 3x3 matrix to 4x4 matrix
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
    d2_t = d2_t @ r_mat
    d2_t = np.array([d2_t[0, 0], d2_t[0, 1], d2_t[0, 2]])  # this is so shit
    d2_t = Dot(coords=d2_t)
    return d2_t


def d_to_ecs(d: Dot, e: Ellipsoid, reverse: bool = None):
    c = e.center
    d_c_t = d.coords - c.coords
    d_t = Dot(coords=d_c_t)
    d_t_r = rot_dot(d_t, e.angle, reverse)
    return d_t_r


def d_in_ell(d: Dot, e: Ellipsoid):
    d = d_to_ecs(d, e)
    x, y, z = d.x, d.y, d.z
    t = e.k * (x ** 2) + e.l * (y ** 2) + e.m * (z ** 2) - e.r ** 2
    return t


def ell_intersection(e1: Ellipsoid, e2: Ellipsoid, precision: float):
    d1 = e1.center
    d2 = e2.center
    if d_in_ell(d1, e2) < precision:
        print('a')
        return True
    d1_t = Dot(0, 0, 0)
    e_temp = Ellipsoid(d1_t, e1.k, e1.l, e1.m, e1.r, e1.angle, -1)

    d2_t = Dot(coords=np.array(d2.coords - d1.coords))

    d2_t = rot_dot(d2_t, e1.angle)

    line = create_line(d1_t, d2_t)
    intersections = line_ellipsoid_intersection(e_temp, line)
    # print(len(intersections), 'intersections')
    for i in intersections:
        i_t = rot_dot(i, e1.angle, True)
        i_t = i_t.coords + e1.center.coords
        # print(i_t)
        # gmsh.model.occ.add_point(i_t[0], i_t[1], i_t[2])
        i_t = Dot(coords=i_t)
        t = d_in_ell(i_t, e2)
        if t < precision:
            print('b')
            return True
    return False



gmsh.initialize()

es = list()

klm = 5
xyz = 2


while len(es) < 23:
    x, y, z = random.random() * klm + 1, random.random() * klm + 1, random.random() * klm + 1
    c = Dot(x, y, z)
    # print('x, y, z =', c.x, ',', c.y, ',', c.z)
    x, y, z = random.random() * np.pi, random.random() * np.pi, random.random() * np.pi
    a = Dot(x, y, z)
    # print('angle =', a.x, ',', a.y, ',', a.z)
    x, y, z = random.random() * xyz + 1, random.random() * xyz + 1, random.random() * xyz + 1
    # print('k, l, m =', x, ',', y, ',', z)
    r = random.random() + 1
    # print('r =', r)
    e = Ellipsoid(c, x, y, z, r, a, len(es))
    f = False
    for e_i in es:
        if ell_intersection(e, e_i, 0.1):
            f = True
    if f:
        continue
    print(f, len(es) + 1)
    e.create_mesh()
    es.append(e)







filename = "GFG.msh"

# origin point
# gmsh.model.occ.add_point(0, 0, 0, 0.001)
#
# x = 2 ** (1/3)
# kek = Ellipsoid(Dot(1, 1, 1), 1, 1, 1, 2, Dot(0, 0, 0), 1)
# kek.create_mesh()
#
# lol = Ellipsoid(Dot(-1, (2 + np.sqrt(2)), 0), 2, 1, 1, 2, Dot(0, np.pi/4, 0), 2)
# lol.create_mesh()
#
# res = ell_intersection(kek, lol)
# print(res)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write(filename)

if 'close' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

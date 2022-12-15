import numpy as np
import geom
import gmsh
import sys


class Line:
    def __init__(self, p: geom.Dot, d: geom.Dot):
        self.p = p
        self.d = d

    def get_dot(self, t: float):
        x = self.p.coords + self.d.coords * t
        x = geom.Dot(coords=x)
        return x


class Ellipsoid:
    def __init__(self, c: geom.Dot, r: float, k: float, l: float, m: float, n: int):
        self.c = c
        self.r = r
        self.k = k
        self.l = l
        self.m = m
        self.n = n

    def create_mesh(self):
        gmsh.model.occ.add_sphere(self.c.x, self.c.y, self.c.z, self.r, tag=self.n)
        gmsh.model.occ.dilate([(3, self.n)], self.c.x, self.c.y, self.c.z, 1/np.sqrt(self.k), 1/np.sqrt(self.l),
                              1/np.sqrt(self.m))


def line_ellipsoid_intersection(e: Ellipsoid, l: Line):
    num_sol = 0
    a = e.k * l.d.x ** 2 + e.l * l.d.y ** 2 + e.m * l.d.z ** 2
    b = 2 * e.k * (l.p.x - e.c.x) * l.d.x + 2 * e.l * (l.p.y - e.c.y) * l.d.y + 2 * e.m * (l.p.z - e.c.z) * l.d.z
    c = e.k * (l.p.x - e.c.x) ** 2 + e.l * (l.p.y - e.c.y) ** 2 + e.m * (l.p.z - e.c.z) ** 2

    if a == 0:
        t = -c/b
        return t
    print(a, b, c)
    dscrm = b ** 2 - a * c * 4

    t = list()

    if dscrm > 0:
        t.append((-b + np.sqrt(dscrm)) / 2 / a)
        t.append((-b + np.sqrt(dscrm)) / 2 / a)
    elif dscrm == 0:
        t.append(-b / 2 / a)

    return t


def create_line(d1: geom.Dot, d2: geom.Dot):
    v = d2.coords - d1.coords
    d = d1
    return Line(d1, geom.Dot(coords=v))


gmsh.initialize()

filename = "GFG.msh"

# origin point
gmsh.model.occ.add_point(0, 0, 0)
gmsh.model.occ.add_point(0, 2, 0)
gmsh.model.occ.add_point(0, (2 + np.sqrt(2)), 0)

kek = Ellipsoid(geom.Dot(0, 0, 0), 2, 1, 1, 2, 1)
kek.create_mesh()

lol = Ellipsoid(geom.Dot(0, (2 + np.sqrt(2)), 0), 2, 2, 2, 1, 2)
lol.create_mesh()

lin = Line(kek.c, geom.Dot(coords=(-kek.c.coords + lol.c.coords)))
print(lin.d.coords)
k = line_ellipsoid_intersection(kek, lin)
print(k)

x = lin.get_dot(k[0])
print(x.coords)

# from Gmsh model.
gmsh.model.geo.synchronize()
gmsh.model.occ.synchronize()

# Generate mesh:
gmsh.model.mesh.generate()

# Write mesh data:
gmsh.write(filename)

# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()

# It finalizes the Gmsh API
gmsh.finalize()
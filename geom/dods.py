import gmsh

from geom.core import *
import numpy as np
import numba
from numba import njit
import time


class Dodecahedron:
    def __init__(self, center: Dot, r: float, angle: Dot, num: int):
        self.center = center
        self.volume = None
        self.surface_loop = None
        self.volume_val = None
        self.r = r
        self.angle = angle
        self.num = num
        self.vertices = []
        self.points = []
        self.pentagons = []
        self.lines = [[0, 0]]
        self.queue = []
        self.normals = []
        self.a = None
        self.r_m = None
        self.build_vertices()

    def build_mesh(self):
        self.build_points()
        self.build_lines()
        self.build_pentagons()

    def build_vertices(self):
        phi = 1 + np.sqrt(5)
        phi = phi / 2
        t = np.sqrt(3)
        self.a = self.r / t / phi * 2
        self.volume_val = np.power(self.a, 3) * (15 + 7 * np.sqrt(5)) / 4
        self.r_m = self.a * phi * phi / 2

        for i in range(0, 2):
            for j in range(0, 2):
                for k in range(0, 2):
                    temp = Dot((1 - 2 * i) * self.r / t, (1 - 2 * j) * self.r / t, (1 - 2 * k) * self.r / t)
                    temp.to_gcs(self.center, self.angle)
                    self.vertices.append(temp)

        for i in range(0, 2):
            for j in range(0, 2):
                temp = Dot(0, (1 - 2 * i) * phi * self.r / t, (1 - 2 * j) * (1 / phi) * self.r / t)
                temp.to_gcs(self.center, self.angle)
                self.vertices.append(temp)

        for i in range(0, 2):
            for j in range(0, 2):
                temp = Dot((1 - 2 * i) * (1 / phi) * self.r / t, 0, (1 - 2 * j) * phi * self.r / t)
                temp.to_gcs(self.center, self.angle)
                self.vertices.append(temp)

        for i in range(0, 2):
            for j in range(0, 2):
                temp = Dot((1 - 2 * i) * phi * self.r / t, (1 - 2 * j) * (1 / phi) * self.r / t, 0)
                temp.to_gcs(self.center, self.angle)
                self.vertices.append(temp)

    def build_points(self):
        lc = 1e-2
        for i in self.vertices:
            temp_point =\
                gmsh.model.occ.add_point(i.x, i.y, i.z)
            self.points.append(temp_point)

    def build_lines(self):
        queue = [[5, 19, 6, 10, 9, 5, 15, 13, 1, 9], [15, 7, 20, 19], [10, 2, 17, 1], [6, 16, 14, 2], [20, 8, 16],
                 [13, 3, 11, 7], [3, 18, 17], [18, 4, 14], [11, 12, 4], [12, 8]]
        t = 0
        for q in queue:
            for i in range(0,
                           len(q) - 1):
                gmsh.model.occ.add_line(q[i] + self.num * 20, q[i + 1] + self.num * 20)
                temp_line = [q[i], q[i + 1]]
                self.lines.append(temp_line)
                t = t + 1

    def build_pentagons(self):
        k = 30 * self.num

        loops = [[1, 2, 3, 4, 5], [5, 6, 7, 8, 9], [9, -4, 13, 14, 15], [-13, -3, 16, 17, 18], [16, -20, -19, 12, 2],
                 [11, 12, -1, 6, 10], [21, -8, 15, 25, 24], [27, 18, 14, -25, 26], [7, 21, 22, 23, - 10],
                 [22, 28, 29, -26, -24], [23, 11, 19, -30, -28], [30, 20, 17, -27, -29]]

        for loop_i in loops:
            line_loops = []
            surfaces = []

            # dots = []
            # dots.append(self.vertices[self.lines[loop_i[0]][0] - 1])

            for line_i in loop_i:
                # dots.append(self.vertices[self.lines[line_i][1] - 1])
                temp_j = line_i
                if temp_j < 0:
                    temp_j = temp_j - k
                else:
                    temp_j = temp_j + k
                line_loops.append(temp_j)

            # self.pentagons.append(inters.Pentagon(dots))

            line1 = self.lines[loop_i[0]]
            line2 = self.lines[loop_i[1]]

            v1 = dist(self.vertices[line1[0] - 1], self.vertices[line1[1] - 1], True)
            v2 = dist(self.vertices[line2[0] - 1], self.vertices[line2[1] - 1], True)

            normal = np.cross(v2, v1)
            normal = normal / np.linalg.norm(normal)
            normal = np.array(normal)
            normal = normal.squeeze()
            normal = Dot(coords=normal)
            self.normals.append(normal)

            curve_loop =\
                gmsh.model.occ.add_curve_loop(line_loops)
            surface = gmsh.model.occ.add_plane_surface([curve_loop])

        surface_loop =\
            gmsh.model.occ.add_surface_loop(range(1 + 12 * self.num, 13 + 12 * self.num))
        self.surface_loop = surface_loop

        # volume = gmsh.model.occ.add_volume([surface_loop])
        # self.volume = volume

    def is_overlapping_with(self, dod: 'Dodecahedron'):
        if dist(dod.center, self.center) <= (dod.r_m + self.r_m):
            return True
        elif (self.r + dod.r) < dist(dod.center, self.center):
            return False
        else:
            for d_i in self.vertices:
                if dist(d_i, dod.center) <= dod.r_m:
                    return True
            for d_i in dod.vertices:
                if dist(d_i, self.center) <= self.r_m:
                    return True

        for n in self.normals:
            mini1, maxi1 = w_interval(self, n)
            mini2, maxi2 = w_interval(dod, n)
            if maxi2 < mini1 or maxi1 < mini2:
                return False

        for n in dod.normals:
            mini1, maxi1 = w_interval(self, n)
            mini2, maxi2 = w_interval(dod, n)
            if maxi2 < mini1 or maxi1 < mini2:
                return False

        for i in self.lines:
            for j in dod.lines:
                line1 = i
                line2 = j
                v1 = dist(self.vertices[line1[0] - 1], self.vertices[line1[1] - 1], True)
                v2 = dist(self.vertices[line2[0] - 1], self.vertices[line2[1] - 1], True)
                n = np.cross(v1, v2)
                n = Dot(coords=n)
                mini1, maxi1 = w_interval(self, n)
                mini2, maxi2 = w_interval(dod, n)
                if maxi2 < mini1 or maxi1 < mini2:
                    return False
        return True

    def in_box(self, d1: Dot, d2: Dot):
        for i in range(3):
            if self.center.coords[i] <= d1.coords[i] or self.center.coords[i] >= d2.coords[i]:
                return False

        for vert in self.vertices:
            for i in range(3):
                if vert.coords[i] <= d1.coords[i] or vert.coords[i] >= d2.coords[i]:
                    return False

        return True


def copy_dod(dod: Dodecahedron, xyz: Dot):
    new_dod = Dodecahedron(xyz, dod.r, dod.angle, dod.num + 1)
    return new_dod


def w_interval(dod: Dodecahedron, n: Dot):
    temp_list = []
    for i in dod.vertices:
        temp_list.append(i.coords)
    dod_coords = np.array([temp_list])
    for i in dod_coords:
        i = np.squeeze(i)
    n_coords = n.coords
    n_coords = np.squeeze(n_coords)
    dod_coords = np.squeeze(dod_coords)
    mini, maxi = interval(dod_coords, n_coords)
    return mini, maxi


@njit
def interval(dod_coords: np.ndarray, n: np.ndarray):
    mini = np.dot(n, dod_coords[0])
    maxi = mini
    for i in dod_coords:
        val = np.dot(n, i)
        if val < mini:
            mini = val
        else:
            maxi = val
    return mini, maxi


def gen_first_dod(first: Dot, second: Dot, angle: Dot, r: float, max_attempts: int = None,
                  epsilon: float = None):
    if max_attempts == 0:
        print("exceeded maximum amount of attempts to build first dod")
        return
    if epsilon is None:
        epsilon = 0

    phi = 1 + np.sqrt(5)
    phi = phi / 2
    t = np.sqrt(3)
    a = r / t / phi * 2
    r_m = a * phi * phi / 2
    x = rand_f(first.x + r + epsilon, second.x + r - epsilon)
    y = rand_f(first.y + r + epsilon, second.y + r - epsilon)
    z = rand_f(first.z + r + epsilon, second.z + r - epsilon)
    center = Dot(x, y, z)
    first_dod = Dodecahedron(center, r, angle, 0)
    if first_dod.in_box(first, second):
        return first_dod
    else:
        return gen_first_dod(first, second, angle, r, max_attempts)


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


def gen_first_dod(first: Dot, second: Dot, angle: Dot, r: float, max_attempts: int = None,
                  epsilon: float = None):
    if max_attempts == 0:
        print("exceeded maximum amount of attempts to build first dod")
        return
    if epsilon is None:
        epsilon = 0

    phi = 1 + np.sqrt(5)
    phi = phi / 2
    t = np.sqrt(3)
    a = r / t / phi * 2
    r_m = a * phi * phi / 2
    x = rand_f(first.x + r + epsilon, second.x + r - epsilon)
    y = rand_f(first.y + r + epsilon, second.y + r - epsilon)
    z = rand_f(first.z + r + epsilon, second.z + r - epsilon)
    center = Dot(x, y, z)
    first_dod = Dodecahedron(center, r, angle, 0)
    if first_dod.in_box(first, second):
        return first_dod
    else:
        return gen_first_dod(first, second, angle, r, max_attempts)


def gen_dod_box(porosity: float, first: Dot, second: Dot, r: float = None,
            r_min: float = None, r_max: float = None, max_time: float = None, max_time_per_1: float = None,
            max_attempts: int = None, epsilon: float = None):
    # first and second dots define box
    # r is radius of circumscribed sphere of dods, it can be None for random value
    # r_min and r_max are borders for random radius
    # max_time is optional limit for total time
    # max_time_per_1 is optional limit for time spent per each dod
    # both are in seconds
    # max_attempts is optional limit for attempts per one dod
    is_random = False

    if r is None:
        is_random = True
    if is_random:
        r = rand_f(r_min, r_max)

    start = time.time()
    target_val = (second.x - first.x) * (second.y - first.y) * (second.z - first.z) * porosity
    attempts = 0

    temp_angle = Dot(is_angle=True)
    dod0 = gen_first_dod(first, second, temp_angle, r, epsilon=epsilon)
    dod0.build_mesh()
    total_volume_val = dod0.volume_val
    dods = list()
    dods.append(dod0)
    surfaces = list()
    surfaces.append(dod0.surface_loop)

    n = 1
    previous = 0

    while attempts < max_attempts and total_volume_val < target_val:
        temp_dot = Dot(border=second)
        temp_angle = Dot(is_angle=True)

        if is_random:
            r = rand_f(r_min, r_max)

        temp_dod = Dodecahedron(temp_dot, r, temp_angle, n)

        intersects = False

        if temp_dod.in_box(first, second) is False:
            intersects = True

        if not intersects:
            for d in dods:
                if temp_dod.is_overlapping_with(d):
                    intersects = True
                    break

        if intersects:
            if max_attempts is not None:
                attempts = attempts + 1
            continue

        attempts = 0

        temp_dod.build_mesh()
        surfaces.append(temp_dod.surface_loop)

        total_volume_val = total_volume_val + temp_dod.volume_val

        dods.append(temp_dod)
        n = n + 1
        print(str(n) + 'th dod made in:')
        now = time.time()
        print('t=' + str(now - start) + 's')

        if max_time_per_1 is not None:
            if max_time_per_1 < now - previous:
                print("time per one limit")
                break

        previous = now

        if max_time is not None:
            if max_time < now - start:
                print("total time limit")
                break

    box_loop = build_box_loop(n, first, second)

    tags = [box_loop]
    tags.extend(surfaces)

    box = gmsh.model.occ.add_volume(tags)
    end = time.time()

if __name__ == '__main__':
    gmsh.initialize()
    gen_dod_box(0.4, Dot(0, 0, 0), Dot(10, 10, 10), r_min=0.5, r_max=2, max_time_per_1=60, max_attempts=100000)
    gmsh.finalize()

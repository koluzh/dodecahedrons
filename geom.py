import numpy as np
from numba import njit
import random
import gmsh
import scipy
# scipy is imported in order to check if its installed because njitted intersection func needs it
# so make sure that you have scipy installed


class Dot:
    def __init__(self, x: float = None, y: float = None, z: float = None, coords: np.array = None,
                 is_angle: bool = None,
                 border: "Dot" = None):

        self.border = border

        if is_angle is None:
            self.is_angle = False
        else:
            self.is_angle = is_angle

        try:
            if coords is not None:
                self.coords = coords
                self.x = self.coords[0, 0]
                self.y = self.coords[0, 1]
                self.z = self.coords[0, 2]

            elif x is not None and y is not None and z is not None:
                self.x = x
                self.y = y
                self.z = z
                self.coords = np.array([x, y, z])
            else:
                self.gen_dot()
        except ValueError:
            print("incorrect arguments")

    def gen_dot(self):
        if self.is_angle:
            self.x = rand_f(0, np.pi * 2)
            self.y = rand_f(0, np.pi * 2)
            self.z = rand_f(0, np.pi * 2)
            self.coords = [self.x, self.y, self.z]
        else:
            self.x = rand_f(0, self.border.x)
            self.y = rand_f(0, self.border.y)
            self.z = rand_f(0, self.border.z)
            self.coords = [self.x, self.y, self.z]

    # to be njitted
    def to_gcs(self, center: 'Dot', angle: 'Dot'):
        t = [[1, 0, 0, center.x],
             [0, 1, 0, center.y],
             [0, 0, 1, center.z],
             [0, 0, 0, 1]]
        t = np.matrix(t, float)

        r1 = [[np.cos(angle.z), np.sin(angle.z), 0],
              [-np.sin(angle.z), np.cos(angle.z), 0],
              [0, 0, 1]]
        r1 = np.array(r1)

        r2 = [[np.cos(angle.y), 0, np.sin(angle.y)],
              [0, 1, 0],
              [-np.sin(angle.y), 0, np.cos(angle.y)]]
        r2 = np.array(r2)

        r3 = [[1, 0, 0],
              [0, np.cos(angle.x), np.sin(angle.x)],
              [0, -np.sin(angle.x), np.cos(angle.x)]]
        r3 = np.array(r3)

        r = np.dot(r1, r2)
        r = np.dot(r, r3)

        self.coords = self.coords @ r

        self.coords = np.array([self.coords[0], self.coords[1], self.coords[2], 1])
        self.coords = t @ self.coords

        self.x = self.coords[0, 0]
        self.y = self.coords[0, 1]
        self.z = self.coords[0, 2]
        self.coords = np.array([self.x, self.y, self.z])


# basically useless
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


def rand_f(start: float, stop: float):
    return random.uniform(start, stop)


def dist(d1: Dot, d2: Dot, is_v: bool = None):
    v = np.matrix(d1.coords) - np.matrix(d2.coords)
    if is_v is None or not is_v:
        length = np.linalg.norm(v)
        return length
    elif is_v:
        return v


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
            temp_point = gmsh.model.geo.add_point(i.x, i.y, i.z)
            self.points.append(temp_point)

    def build_lines(self):
        queue = [[5, 19, 6, 10, 9, 5, 15, 13, 1, 9], [15, 7, 20, 19], [10, 2, 17, 1], [6, 16, 14, 2], [20, 8, 16],
                 [13, 3, 11, 7], [3, 18, 17], [18, 4, 14], [11, 12, 4], [12, 8]]
        t = 0
        for q in queue:
            for i in range(0, len(q) - 1):
                gmsh.model.geo.add_line(q[i] + self.num * 20, q[i + 1] + self.num * 20)
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
            normal = Dot(coords=normal)
            self.normals.append(normal)

            curve_loop = gmsh.model.geo.add_curve_loop(line_loops)
            surface = gmsh.model.geo.add_plane_surface([curve_loop])

        surface_loop = gmsh.model.geo.add_surface_loop(range(1 + 12 * self.num, 13 + 12 * self.num))
        self.surface_loop = surface_loop
        # volume = gmsh.model.geo.add_volume([surface_loop])
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


@njit()
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


def build_box_loop(num: int, s: Dot, d: Dot):
    line_n = 30 * num
    point_n = 20 * num
    surface_n = 12 * num
    for i in range(2):
        for j in range(2):
            for k in range(2):
                gmsh.model.geo.add_point(s.coords[0] + d.coords[0] * i, s.coords[1] + d.coords[1] * j,
                                         s.coords[2] + d.coords[2] * k)
    queue = [[1, 2, 6, 5, 1], [2, 4, 3, 1], [6, 8, 4], [5, 7, 8], [7, 3]]
    lines = list()

    for q in queue:
        for i in range(0, len(q) - 1):
            gmsh.model.geo.add_line(q[i] + point_n, q[i + 1] + point_n)
            temp_line = [q[i], q[i + 1]]
            lines.append(temp_line)

    line_queue = [[1, 2, 3, 4], [8, 9, -5, 2], [10, 11, -8, 3], [6, -12, 11, 9], [10, 12, 7, -4], [1, 5, 6, 7]]

    for loop_i in line_queue:
        line_loops = []

        for line_i in loop_i:
            # dots.append(self.vertices[self.lines[line_i][1] - 1])
            temp_j = line_i
            if temp_j < 0:
                temp_j = temp_j - line_n
            else:
                temp_j = temp_j + line_n
            line_loops.append(temp_j)

        curve_loop = gmsh.model.geo.add_curve_loop(line_loops)
        surface = gmsh.model.geo.add_plane_surface([curve_loop])

    surface_loop = gmsh.model.geo.add_surface_loop(range(1 + surface_n, 7 + surface_n))
    return surface_loop

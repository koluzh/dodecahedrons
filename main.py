import numpy as np
import gmsh
import sys
import random
from numba import njit
import time
# you also need to install scipy 0.16+
# import intersection as inters
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

if __name__ == "__main__":
    gmsh.initialize()

VOLUMES = []


def rand_f(start: float, stop: float):
    return random.uniform(start, stop)


class Dot:
    def __init__(self, x: float = None, y: float = None, z: float = None, coords: np.array = None, is_angle: bool = None,
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
        phi = phi/2
        t = np.sqrt(3)
        self.a = self.r / t / phi * 2
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
            for i in range(0, len(q)-1):
                gmsh.model.geo.add_line(q[i] + self.num * 20, q[i + 1] + self.num * 20)
                temp_line = [q[i], q[i+1]]
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
            normal = normal/np.linalg.norm(normal)
            normal = Dot(coords=normal)
            self.normals.append(normal)

            curve_loop = gmsh.model.geo.add_curve_loop(line_loops)
            surface = gmsh.model.geo.add_plane_surface([curve_loop])

        surface_loop = gmsh.model.geo.add_surface_loop(range(1 + 12 * self.num, 13 + 12 * self.num))
        volume = gmsh.model.geo.add_volume([surface_loop])
        VOLUMES.append(volume)

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


if __name__ == '__main__':
    theta = Dot(0, 0, 0)
    BORDER = Dot(10, 10, 10)
    temp_dot = Dot(border=BORDER)
    temp_angle = Dot(is_angle=True)
    MAX_ATTEMPTS = 10000

    start = time.time()
    N = 1
    dods = []
    dod0 = Dodecahedron(theta, 2, theta, 0)
    dod0.build_mesh()
    dods.append(dod0)
    attempts = 0
    # parallelize
    while attempts < MAX_ATTEMPTS:
        temp_dot = Dot(border=BORDER)
        temp_angle = Dot(is_angle=True)

        # print(temp_dot.coords)
        # print(temp_angle.coords)
        # print(N)

        temp_dod = Dodecahedron(temp_dot, 2, temp_angle, N)

        intersects = False

        for d in dods:
            if temp_dod.is_overlapping_with(d):
                intersects = True
                break

        if intersects:
            attempts = attempts + 1
            # print("skipped")
            continue

        attempts = 0

        temp_dod.build_mesh()

        dods.append(temp_dod)
        N = N + 1

    end = time.time()
    print(end - start)
    print(len(VOLUMES))
    # Create the relevant Gmsh data structures
    # from Gmsh model.
    gmsh.model.geo.synchronize()

    # Generate mesh:
    gmsh.model.mesh.generate()

    # Write mesh data:
    gmsh.write("GFG.msh")

    # Creates  graphical user interface
    if 'close' not in sys.argv:
        gmsh.fltk.run()

    # It finalizes the Gmsh API
    gmsh.finalize()
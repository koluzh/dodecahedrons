import numpy as np
import gmsh
import random
from scipy.spatial.transform import Rotation as rot
from typing import Union

class Dot:
    def __init__(self, x: float = None, y: float = None, z: float = None,
                 coords: np.ndarray = None):
        if not any(i is None for i in [x, y, z]):
            self._x = x
            self._y = y
            self._z = z
            self._coords = np.array([x, y, z])
        elif coords is not None:
            self.check_vector(coords)
            self._coords = np.squeeze(np.asarray(coords))[:3]
            self._x = coords[0]
            self._y = coords[1]
            self._z = coords[2]

    def check_vector(self, coords: np.ndarray) -> None:
        if not isinstance(coords, np.ndarray):
            raise Exception("not numpy array")
        if coords.size < 3 or coords.size > 4:
            raise Exception("incorrect size of vector")

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @x.setter
    def x(self, value: float):
        self._x = value
        self.coords[0] = value

    @y.setter
    def y(self, value: float):
        self._y = value
        self.coords[1] = value

    @z.setter
    def z(self, value: float):
        self._z = value
        self.coords[2] = value

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, value: np.array):
        self.check_vector(value)
        self._coords = np.squeeze(np.asarray(value))[:3]

    # to be njitted
    def to_gcs(self, center: 'Dot', angle: 'Dot'):
        # needed to convert coordinates to gcs from local coordinate system
        t = [[1, 0, 0, center.x],
             [0, 1, 0, center.y],
             [0, 0, 1, center.z],
             [0, 0, 0, 1]]
        t = np.matrix(t, float)

        r1 = [[np.cos(angle.z),  np.sin(angle.z),  0],
              [-np.sin(angle.z), np.cos(angle.z), 0],
              [0,                0,                1]]
        r1 = np.array(r1)

        r2 = [[np.cos(angle.y),  0, np.sin(angle.y)],
              [0,                1, 0],
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
        self.coords = t @ np.append(self.coords, [1])

        self.x = self.coords[0]
        self.y = self.coords[1]
        self.z = self.coords[2]
        self.coords = np.array([self.x, self.y, self.z])

    def as_list(self):
        coord_list = [self.coords[0], self.coords[1], self.coords[2]]
        return coord_list

    def build_mesh(self):
        gmsh.model.occ.add_point(self.x, self.y, self.z)


def rand_f(start: float, stop: float):
    return random.uniform(start, stop)


def dist(d1: Dot, d2: Dot, is_v: bool = None):
    v = np.matrix(d1.coords) - np.matrix(d2.coords)
    if not is_v:
        length = np.linalg.norm(v)
        return length
    elif is_v:
        return v


def get_rotation_matrix(angle: Dot):
    # print(angle.coords)
    r = rot.from_euler('xyz', angle.coords)
    # print(r.as_matrix())
    return r.as_matrix()


def build_box_loop(num: int, s: Dot, d: Dot):
    line_n = 30 * num
    point_n = 20 * num
    surface_n = 12 * num
    for i in range(2):
        for j in range(2):
            for k in range(2):
                gmsh.model.occ.add_point(s.coords[0] + d.coords[0] * i, s.coords[1] + d.coords[1] * j,
                                         s.coords[2] + d.coords[2] * k)
    queue = [[1, 2, 6, 5, 1], [2, 4, 3, 1], [6, 8, 4], [5, 7, 8], [7, 3]]
    lines = list()

    for q in queue:
        for i in range(0, len(q) - 1):
            gmsh.model.occ.add_line(q[i] + point_n, q[i + 1] + point_n)
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

        curve_loop =\
            gmsh.model.\
                occ.add_curve_loop(line_loops)
        surface = gmsh.model.occ.add_plane_surface([curve_loop])

    surface_loop =\
        gmsh.model.occ.add_surface_loop(range(1 + surface_n, 7 + surface_n))
    return surface_loop


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


class Box:
    def __init__(self, s: Dot, d: Dot):
        self.s = s
        self.d = d
        self.c = Dot(coords=(s.coords + d.coords) / 2)   # center
        dots = list()
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    temp_dot = Dot(s.coords[0] + (d.coords[0] - s.coords[0]) * i,
                                             s.coords[1] + (d.coords[1] - s.coords[1]) * j,
                                             s.coords[2] + (d.coords[2] - s.coords[2]) * k)
                    dots.append(temp_dot)
                    gmsh.model.occ.add_point(s.coords[0] + (d.coords[0] - s.coords[0]) * i,
                                             s.coords[1] + (d.coords[1] - s.coords[1]) * j,
                                             s.coords[2] + (d.coords[2] - s.coords[2]) * k)
        queue = [[1, 2, 6, 5, 1], [2, 4, 3, 1], [6, 8, 4], [5, 7, 8], [7, 3]]   # queue of dots to connect into lines
        plane_nums = [[1, 2, 6, 5], [2, 4, 3, 1], [2, 4, 8, 6], [6, 8, 7, 5], [5, 7, 3, 1], [8, 4, 3, 7]]
        planes = list()
        # x: 5786, 1342 plane_nums[3], plane_nums[1]
        # y: 4378, 2156 plane_nums[0], planes_nums[5]
        # z: plane_nums[2], plane_nums[4]
        for nums in plane_nums:
            plane_dots = list()
            for i in nums:
                plane_dots.append(dots[i - 1])
            temp_plane = Plane(plane_dots, self.c)
            planes.append(temp_plane)

        self.planes = planes

        lines = list()

        for q in queue:
            for i in range(0, len(q) - 1):
                gmsh.model.occ.add_line(q[i], q[i + 1])
                temp_line = [q[i], q[i + 1]]
                lines.append(temp_line)

    def build_mesh(self, n):
        line_queue = [[1, 2, 3, 4], [8, 9, -5, 2], [10, 11, -8, 3], [6, -12, 11, 9], [10, 12, 7, -4], [1, 5, 6, 7]]

        for loop_i in line_queue:
            line_loops = []

            for line_i in loop_i:
                # dots.append(self.vertices[self.lines[line_i][1] - 1])
                temp_j = line_i
                if temp_j < 0:
                    temp_j = temp_j
                else:
                    temp_j = temp_j
                line_loops.append(temp_j)

            curve_loop = gmsh.model.occ.add_curve_loop(line_loops)
            surface = gmsh.model.occ.add_plane_surface([curve_loop])

        surface_loop = gmsh.model.occ.add_surface_loop(range(n + 1, n + 7))
        return surface_loop


class Plane:
    def __init__(self, d: list[Dot], c: Dot = None):
        # d is list of dots to define normal vector of a plane and its D
        # c is center of
        if len(d) < 3:
            print('Dot amount in plane is wrong')
        self.dots_list = d
        v1 = np.array(d[0].coords - d[1].coords)
        v2 = np.array(d[2].coords - d[1].coords)
        self.normal = np.cross(v1, v2)
        self.normal = self.normal / np.linalg.norm(self.normal)
        self.d = np.dot(self.normal, self.dots_list[0].coords)
        if c is not None:
            self.sign = np.dot(self.normal, c.coords)
        else:
            self.sign = 0

    def get_projection(self, dot: Dot):
        # print(dot.coords)
        d = np.dot(dot.coords - self.dots_list[0].coords, self.normal)
        if d > 0:
            q_c = dot.coords - d * self.normal
        else:
            q_c = dot.coords + d * self.normal
        q = Dot(coords=q_c)
        # print(q.coords)
        return q

if __name__ == '__main__':
    d1 = Dot(1, 2, 3)
    d2 = Dot(2, 2, 2)

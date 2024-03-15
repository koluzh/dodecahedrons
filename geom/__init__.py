from . import dods
from . import core
from .core import *
from . import ellps
from . import sphrs
from . import cyls
from typing import Union, Type
import time
import gmsh


class PorousBox:
    class State:
        def __init__(self, pb: 'PorousBox'):
            self.starting_time = time.time()
            self.cur_obj_attempts = 0
            self.cur_obj_start_time = time.time()
            self.num_of_objs = 0
            self.total_attempts = 0
            self.__parent = pb
            self.cur_porocity = 0
            self.cur_void_volume = 0

        def attempt(self):
            self.cur_obj_attempts = 0
            self.cur_obj_start_time = time.time()

        def success(self, obj):
            self.cur_void_volume += obj.volume
            self.cur_porocity = self.cur_void_volume / self.__parent.box.volume

        @property
        def total_time(self):
            return time.time() - self.starting_time

        @property
        def cur_obj_time(self):
            return time.time() - self.cur_obj_start_time


    def __init__(self,
                 box_a: Dot,
                 box_b: Dot,
                 target_porosity: float,
                 obj_type: Union[Type[ellps.Ellipsoid] | Type[dods.Dodecahedron]],
                 r_min: float, r_max: float,
                 stop_func: callable = None,
                 **kwargs):
        '''
        :param box_a:
            s param of Box class
        :param box_b:
            d param of Box class
        :param target_porosity:
            target porosity to be achieved by this class
        :param obj_type:
            Ellipsoid or Dodecahedron, later Cylinder and Sphere
        :param r_min:
            minimum outer radius of objects
        :param r_max:
            maximum outer radius of objects
        :param stop_func:

        :param kwargs:
            ells:
                k_min, k_max
                l_min, l_max
                m_min, m_max
            dods:
                no specific kwargs

        '''
        self.box = Box(box_a, box_b)
        self.obj_type = obj_type
        self.r_min = r_min
        self.r_max = r_max
        self.kwargs = kwargs
        self.cur_num = 0
        self.target_porosity = target_porosity
        self.state = self.State(self)
        self.stop = lambda state: state.cur_porocity >= self.target_porosity if stop_func is None else stop_func
        self.objects = []
        self.s = box_a
        self.d = box_b

    def get_rand_dot_in_box(self, r: float = None):
        if r is None:
            r = 0
        start = self.box.s
        stop = self.box.d
        x = rand_f(start.x + r, stop.x - r)
        y = rand_f(start.y + r, stop.y - r)
        z = rand_f(start.z + r, stop.z - r)
        return Dot(x, y, z)

    def rand_r(self) -> float:
        return rand_f(self.r_min, self.r_max)

    def rand_angle(self) -> Dot:
        angle_x = rand_f(-np.pi / 2, np.pi / 2)
        angle_y = rand_f(-np.pi / 2, np.pi / 2)
        angle_z = rand_f(-np.pi / 2, np.pi / 2)
        return Dot(angle_x, angle_y, angle_z)

    def gen_rand_ell(self) -> ellps.Ellipsoid:
        k_min = self.kwargs['k_min']
        k_max = self.kwargs['k_max']
        l_min = self.kwargs['l_min']
        l_max = self.kwargs['l_max']
        m_min = self.kwargs['m_min']
        m_max = self.kwargs['m_max']
        k = rand_f(k_min, k_max)
        l = rand_f(l_min, l_max)
        m = rand_f(m_min, m_max)
        r = self.rand_r()
        a = r / np.sqrt(k)
        b = r / np.sqrt(l)
        c = r / np.sqrt(m)
        center = self.get_rand_dot_in_box(
            r=max([a, b, c])
        )
        angle = self.rand_angle()
        self.cur_num += 1
        return ellps.Ellipsoid(center, k, l, m, r, angle, self.cur_num)

    def gen_rand_dod(self) -> dods.Dodecahedron:
        self.cur_num += 1
        return dods.Dodecahedron(
            self.get_rand_dot_in_box(),
            self.rand_r(),
            self.rand_angle(),
            self.cur_num
        )

    def gen_rand_obj(self):
        result = None
        match self.obj_type:
            case ellps.Ellipsoid:
                result = self.gen_rand_ell()
            case dods.Dodecahedron:
                result = self.gen_rand_dod()
        return result

    def get_obj_in_box(self):
        # self.state.update() - updates time and cur and total attempts
        # self.state.update(dv=some_val) - updates attempts (sets current obj attempts to zero, adds to total),
        # updates void volume, porocity, updates time
        while not self.stop(self.state):
            obj = self.gen_rand_obj()
            if obj.in_box(self.box):
                return obj

    def build_box_loop(self, num: int, s: Dot, d: Dot) -> None:
        """
            :param num:
            ;id number of this 3d object in gmsh
            :param s:
            ;starting dot
            :param d:
            ;dot that represents sides length (i know this is stupid af)
            :return:
        """
        # line_n = 30 * num
        # point_n = 20 * num
        # surface_n = 12 * num
        line_n = 0
        point_n = 0
        surface_n = 0
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
        dots = []
        for loop_i in line_queue:
            line_loops = []

            for line_i in loop_i:
                # dots.append(vertices[lines[line_i][1] - 1])
                temp_j = line_i
                if temp_j < 0:
                    temp_j = temp_j - line_n
                else:
                    temp_j = temp_j + line_n
                line_loops.append(temp_j)

            curve_loop = gmsh.model.occ. \
                add_curve_loop(line_loops)
            surface = gmsh.model.occ.add_plane_surface([curve_loop])

        surface_loop = \
            gmsh.model.occ.add_surface_loop(range(1 + surface_n, 7 + surface_n))
        return surface_loop

    def does_not_intersect_others(self, obj):
        for o in self.objects:
            if obj.intersects(o):
                return False
        return True

    def fill_box(self):
        while not self.stop(self.state):
            # print("not stop")
            self.state.attempt()
            obj = self.get_obj_in_box()
            print("got obj in box")
            if self.does_not_intersect_others(obj):
                print('yay!')
                obj.num = len(self.objects)
                print(self.state.cur_porocity)
                obj.build_mesh()
                self.objects.append(obj)
                self.state.success(obj)
        n = len(self.objects)
        box_volume = gmsh.model.occ.addBox(self.s.x, self.s.y, self.s.z, self.d.x - self.s.x, self.d.y - self.s.y, self.d.z - self.s.z)
        box_loops = gmsh.model.occ.getSurfaceLoops(box_volume)
        print(f'box loops: {box_loops}')
        print(f'box_volume: {box_volume}')
        gmsh.model.occ.remove([(3, x) for x in range(n)])
        gmsh.model.occ.remove([(3, box_volume)])
        tags = []
        tags = list(box_loops[0])
        tags.extend(x + 1 for x in range(n))
        print(tags)
        box = gmsh.model.occ.add_volume(tags)



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
    # print(r)

    r_mat = r
    d2_t = d2_t @ r_mat
    d2_t = np.array([d2_t[0, 0], d2_t[0, 1], d2_t[0, 2]])  # this is so shit
    d2_t = Dot(coords=d2_t)
    return d2_t


def d_to_ecs(d: Dot, e: Ellipsoid, reverse: bool = None):
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


def ell_intersection(e1: Ellipsoid, e2: Ellipsoid, epsilon: float = None, testing: bool = None):
    if epsilon is None:
        epsilon = 0
    # E1 CENTER MUST BE 0,0,0 I.E.
    d1 = e1.center
    d2 = e2.center
    d1_t = Dot(0, 0, 0)
    t = d_in_ell(d1, e2)
    if t > 0:
        t = np.sqrt(t)
    if t <= epsilon:
        return True, list()
    e_temp = Ellipsoid(d1_t, e1.k, e1.l, e1.m, e1.r, e1.angle, -1)
    d2_t = Dot(coords=np.array(np.array(d2.coords) - np.array(d1.coords)))
    d2_t = rot_dot(d2_t, e1.angle)
    line = create_line(d1_t, d2_t)
    intersections = line_ellipsoid_intersection(e_temp, line)
    inter_dots = list()
    for i in intersections:
        i_t = rot_dot(i, e1.angle, reverse=True)
        i_t = i_t.coords + e1.center.coords
        # print(i_t)
        # gmsh.model.occ.add_point(i_t[0], i_t[1], i_t[2])
        i_t = Dot(coords=i_t)
        inter_dots.append(i_t)
        t = d_in_ell(i_t, e2)
        if t > 0:
            t = np.sqrt(t)
        if t <= epsilon:
            return True, inter_dots
    return False, inter_dots


def gen_first_ell(box: Box, angle: Dot, r: float, k: float, l: float, m: float, max_attempts: int = None,
                  epsilon: float = None, n=0):
    if max_attempts == 0:
        print("exceeded maximum amount of attempts to build first ell")
        return
    if epsilon is None:
        epsilon = 0

    first = box.s
    second = box.d

    x = rand_f(first.x + r + epsilon, second.x + r - epsilon)
    y = rand_f(first.y + r + epsilon, second.y + r - epsilon)
    z = rand_f(first.z + r + epsilon, second.z + r - epsilon)

    # print(x, y, z)

    center = Dot(x, y, z)
    first_ell = Ellipsoid(center, k, l, m, r, angle, n)
    ###
    # first_ell.create_mesh()
    ###
    if first_ell.in_box(box, epsilon):
        return first_ell
    else:
        if max_attempts is not None:
            max_attempts = max_attempts - 1
        return gen_first_ell(box, angle, r, k, l, m, max_attempts)


def gen_ell_box(porosity: float, first: Dot, second: Dot, k: float = None,
                l: float = None, m: float = None, r: float = None, k_min: float = None,
                k_max: float = None, r_min: float = None, r_max: float = None, max_time: float = None,
                max_time_per_1: float = None, max_attempts: int = None, epsilon: float = None, testing: bool = None):
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
        k = rand_f(k_min, k_max)
        l = rand_f(k_min, k_max)
        m = rand_f(k_min, k_max)

    start = time.time()
    box_val = (second.x - first.x) * (second.y - first.y) * (second.z - first.z)
    target_val = box_val * porosity
    attempts = 0
    total_attempts = 0

    box = Box(first, second)

    temp_angle = Dot(is_angle=True)
    ell0 = gen_first_ell(box, temp_angle, r, k, l, m, max_attempts, epsilon=epsilon)

    if ell0 is None:
        return

    ell0.create_mesh()
    total_volume_val = ell0.volume
    if testing is True:
        ell0.print_info()

    ells = list()
    ells.append(ell0)

    n = 1
    previous = 0

    inter_dots_big = list()

    while attempts < max_attempts and total_volume_val < target_val:
        temp_dot = Dot(border=second)
        temp_angle = Dot(is_angle=True)

        if is_random:
            r = rand_f(r_min, r_max)
            k = rand_f(k_min, k_max)
            l = rand_f(k_min, k_max)
            m = rand_f(k_min, k_max)

        temp_ell = Ellipsoid(temp_dot, k, l, m, r, temp_angle, n)

        intersects = False

        if temp_ell.in_box(box, epsilon=epsilon) is False:
            intersects = True

        inter_dots = []

        if not intersects:
            for e in ells:
                inter_bool, inter_dots = ell_intersection(temp_ell, e, epsilon=epsilon, testing=testing)
                if inter_bool:
                    intersects = True
                    break

        if intersects:
            if max_attempts is not None:
                attempts = attempts + 1
            continue

        inter_dots_big.extend(inter_dots)

        total_attempts = total_attempts + attempts

        attempts = 0

        temp_ell.create_mesh()

        total_volume_val = total_volume_val + temp_ell.volume

        ells.append(temp_ell)
        if testing is True:
            temp_ell.print_info()
        n = n + 1
        print(str(n) + 'th ell made in:')
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
    if attempts >= max_attempts:
        print('out of attempts\n', 'porosity =', total_volume_val/box_val)
    else:
        print('succesful')
    if testing:
        for d in inter_dots_big:
            if type(d) is not Dot:
                continue
            gmsh.model.occ.add_point(d.x, d.y, d.z)
    print(total_attempts)

    #gmsh.model.occ.cut()
    end = time.time()
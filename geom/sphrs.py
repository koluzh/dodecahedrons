from __init__ import *


# ??? replace dods with spheres

def gen_sphere_box(porosity: float, first: Dot, second: Dot, r: float = None,
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
        print(str(n) + 'th sphere made in:')
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

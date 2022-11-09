import gmsh
import sys
import time
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
import geom

if __name__ == "__main__":
    gmsh.initialize()

VOLUMES = []


if __name__ == '__main__':
    theta = geom.Dot(0, 0, 0)
    BORDER = geom.Dot(10, 15, 20)
    temp_dot = geom.Dot(border=BORDER)
    temp_angle = geom.Dot(is_angle=True)
    MAX_ATTEMPTS = 10000

    start = time.time()

    # box

    # porosity
    porosity = 0.4
    target_val = porosity * 1000

    # first dod characteristics
    dod0 = geom.Dodecahedron(geom.Dot(5, 5, 5), 2, theta, 0)
    dod0.build_mesh()

    # starting settings
    N = 1
    dods = list()
    SURFACES = list()
    dods.append(dod0)
    total_volume_val = dod0.volume_val
    SURFACES.append(dod0.surface_loop)
    # VOLUMES.append(dod0.volume)
    attempts = 0

    while attempts < MAX_ATTEMPTS and total_volume_val < target_val:
        temp_dot = geom.Dot(border=BORDER)
        temp_angle = geom.Dot(is_angle=True)

        # print(temp_dot.coords)
        # print(temp_angle.coords)
        # print(N)

        temp_dod = geom.Dodecahedron(temp_dot, 2, temp_angle, N)

        intersects = False

        if temp_dod.in_box(theta, BORDER) is False:
            intersects = True

        if intersects is False:
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
        # VOLUMES.append(temp_dod.volume)
        SURFACES.append(temp_dod.surface_loop)

        total_volume_val = total_volume_val + temp_dod.volume_val

        dods.append(temp_dod)
        N = N + 1

    box_loop = geom.build_box_loop(N, theta, BORDER)

    tags = [box_loop]
    tags.extend(SURFACES)
    box = gmsh.model.geo.add_volume(tags)

    end = time.time()
    print(end - start)
    print(SURFACES[0])
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

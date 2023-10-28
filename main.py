import gmsh
import sys
import geom
from geom.core import *
import geom.dods
import numpy as np
import time


if __name__ == '__main__':
    gmsh.initialize()

    # settings

    # precision
    EPSILON = 1e-6

    # two dots defining box
    # might not work if theta is not (0,0,0)
    theta = Dot(0, 0, 0)
    BORDER = Dot(10, 15, 20)

    # maximum amount of attempts per 1 dod, after which gen_box will jump to creating box volume, can be None
    maximum = 10000

    # maximum amount of total time in seconds allowed to generate dods,
    # after which gen_box will jump to creating box volume, can be None
    maximum_time = 3600

    # maximum amount of time in seconds allowed to generate ONE dod,
    # after which gen_box will jump to creating box volume, can be None
    max_time_per_one = 60

    # radius of circumscribed sphere of dodecahedrons to be built
    # can be None
    radius = 2

    # in case of radius being None you need to define these
    r_min = None
    r_max = None

    filename = "GFG.msh"
    start = time.time()

    # function that creates box volume with cavities
    #
    # kek = Box(Dot(0, 0, 0), Dot(10, 10, 10))

    geom.dods.gen_dod_box(0.4, Dot(0, 0, 0), Dot(10, 10, 10), 2, max_attempts=100000)

    end = time.time()
    print('time', end - start)

    # from Gmsh model.
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

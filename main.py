import gmsh
import sys
import geom


if __name__ == '__main__':
    gmsh.initialize()

    # settings

    # precision
    EPSILON = 1e-6

    # two dots defining box
    # might not work if theta is not (0,0,0)
    theta = geom.Dot(0, 0, 0)
    BORDER = geom.Dot(10, 15, 20)

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

    porosity = 0.1

    # function that creates box volume with cavities
    # geom.gen_box(porosity, theta, BORDER, r=radius,
    #             r_min=r_min, r_max=r_max, max_time=maximum_time, max_attempts=maximum, epsilon=EPSILON)

    tag = gmsh.model.occ.add_sphere(0, 0, 0, 5, 1)

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

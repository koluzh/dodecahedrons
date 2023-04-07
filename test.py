import numpy as np
from geom import *
import gmsh
import random
import sys

filename = "GFG.msh"
gmsh.initialize()

es = list()
# origin point

dot0 = Dot(0, 0, 0).build_mesh()
dot1 = Dot(1, 0, 0).build_mesh()
dot2 = Dot(0, 1, 0).build_mesh()
dot3 = Dot(0, 0, 1).build_mesh()

from geom.cyls import Cylinder

c1 = Cylinder(Dot(0, 0, 0), 2, 5, 1)
c1.build_mesh()

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write(filename)

if 'close' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

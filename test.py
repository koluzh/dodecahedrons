import numpy as np
from geom import *
import gmsh
import random
import sys





gmsh.initialize()

es = list()


filename = "GFG.msh"

# origin point

dot1 = Dot(0, 0, 0)
dot2 = Dot(10, 10, 10)

box = Box(dot1, dot2)

xyz = [1.2571628924859857, 8.418863384407327, 8.18503094890025]

k, l, m = [1.4575589628135759, 3.8813418966060844, 1.3417454290908415]

angle = [0.3212844891761337, 0.37135220813792363, 1.0243959529595408]
kek = Ellipsoid(Dot(coords=np.array(xyz)), k, l, m, 2,
                Dot(coords=np.array(angle)), 1)
kek.in_box(box, 0.1)

kek.create_mesh()

# xyz = [2.8861338297589296, 1.7990208685564153, 1.3368151487506896]
# k, l, m = [3.9200293893628655, 1.5705003060664666, 3.304796015160285]
# angle = [5.285209712343209, 0.6781614221526329, 0.549643614312186]
# lol = Ellipsoid(Dot(coords=np.array(xyz)), k, l, m, 2,
#                 Dot(coords=np.array(angle)), 2)
# lol.create_mesh()

# res = ell_intersection(kek, lol, testing=True)
# print(res)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write(filename)

if 'close' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

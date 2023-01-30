import numpy as np
import geom
import gmsh
import sys
filename = 'GFG.msh'

dot1 = geom.Dot(0, 0, 0)
dot2 = geom.Dot(10, 10, 10)
por = 0.4



gmsh.initialize()

# box = geom.Box(dot1, dot2)
#
# gmsh.model.occ.add_point(box.c.coords[0], box.c.coords[1], box.c.coords[2])
# q = box.planes[4].get_projection(box.c)
# gmsh.model.occ.add_point(q[0], q[1], q[2])

geom.gen_ell_box(0.5, dot1, dot2, k_min=1, k_max=4, r_min=1, r_max=2, max_attempts=500)


gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write(filename)

if 'close' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

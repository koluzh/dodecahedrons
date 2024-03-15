from geom.stop_condition import get_stop_function
from geom import PorousBox
from geom.core import Dot
from geom.ellps import Ellipsoid
from geom.dods import Dodecahedron
import gmsh
import sys

a = Dot(0, 0, 0)
b = Dot(10, 10, 10)

stop_f = get_stop_function(0.2)

filename = 'GFG.msh'
gmsh.initialize()

pb = PorousBox(
    box_a=a, box_b=b, target_porosity=0.15, obj_type=Dodecahedron, r_min=1, r_max=2,
    k_min=1, k_max=2,
    l_min=1, l_max=2,
    m_min=1, m_max=2
)
pb.fill_box()


gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write(filename)

if 'close' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

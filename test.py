import numpy as np
from geom.ellps import Ellipsoid, d_in_ell, ell_2_ell
from geom.core import get_rotation_matrix, Dot
import gmsh
import sys

from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def ellipsoid_intersection_test(Sigma_A, Sigma_B, mu_A, mu_B, tau):
    eigvals, eigvecs = np.linalg.eig(Sigma_B)
    l = np.diag(eigvals)
    print(l)
    lambdas, Phi, v_squared = ellipsoid_intersection_test_helper(Sigma_A, l, mu_A, mu_B)
    res = minimize_scalar(ellipsoid_K_function,
                          bracket=[0.0, 0.5, 1.0],
                          args=(lambdas, v_squared, tau))
    return (res.fun[0] >= 0 - 0.000_001)


def ellipsoid_intersection_test_helper(Sigma_A, Sigma_B, mu_A, mu_B):
    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    return lambdas, Phi, v_squared


def ellipsoid_K_function(ss, lambdas, v_squared, tau):
    ss = np.array(ss).reshape((-1,1))
    lambdas = np.array(lambdas).reshape((1,-1))
    v_squared = np.array(v_squared).reshape((1,-1))
    return 1.-(1./tau**2)*np.sum(v_squared*((ss*(1.-ss))/(1.+ss*(lambdas-1.))), axis=1)


if __name__ == '__main__':
    gmsh.initialize()
    d1 = Dot(0, 0, 0)
    d2 = Dot(4, 0, 0)
    angle = Dot(np.pi/4, np.pi/8, 0)
    E1 = Ellipsoid(d1, 1, 1, 1, 2, d1, 1)
    # 2 1 1
    E1.create_mesh()
    E2 = Ellipsoid(d2, 1, 1, 1, 2, d1, 2)
    # 3 2 4
    E2.create_mesh()
    print(E2.l)
    d = np.array([2, 0, 0]) - E2.center.coords
    print(d @ E2.matrix @ d.T)
    print(ellipsoid_intersection_test(np.linalg.inv(E1.matrix), np.linalg.inv(E2.matrix), E1.center.coords, E2.center.coords, tau=1))
    # print("eigvals:")
    # for i in eigvals:
    #     print(i)
    # print("\neigvecs:")
    # for i in eigvecs:
    #     print(i)
    # print("\n")
    # multiply vector by some factor so that the 4th component equals 1,
    # check if point represented by this vector is inside of both ellipsoids

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate()
    gmsh.write('GFG.msh')

    if 'close' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

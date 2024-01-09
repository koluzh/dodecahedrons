import numpy as np
from geom.ellps import Ellipsoid, d_in_ell
from geom.core import get_rotation_matrix, Dot
import gmsh
import sys

def Algebraic_Separation_Condition(coeff_canon_i, coeff_canon_j, r_i, r_j, A_i, A_j):
    # Ellipsoid matrices in the canonical form.s
    A = np.array([[1 / coeff_canon_i[0] ** 2, 0, 0, 0],
                  [0, 1 / coeff_canon_i[1] ** 2, 0, 0],
                  [0, 0, 1 / coeff_canon_i[2] ** 2, 0],
                  [0, 0, 0, -1]])
    B = np.array([[1 / coeff_canon_j[0] ** 2, 0, 0, 0],
                  [0, 1 / coeff_canon_j[1] ** 2, 0, 0],
                  [0, 0, 1 / coeff_canon_j[2] ** 2, 0],
                  [0, 0, 0, -1]])

    # Rigid body transformations.
    T_i = np.vstack((np.hstack((A_i, r_i)), np.array([0, 0, 0, 1])))
    T_j = np.vstack((np.hstack((A_j, r_j)), np.array([0, 0, 0, 1])))
    Ma = T_i
    Mb = T_j

    # aij belongs to A in det(lambda*A - Ma'*(Mb^-1)'*B*(Mb^-1)*Ma).
    a11, a12, a13, a14 = A[0, 0], A[0, 1], A[0, 2], A[0, 3]
    a21, a22, a23, a24 = A[1, 0], A[1, 1], A[1, 2], A[1, 3]
    a31, a32, a33, a34 = A[2, 0], A[2, 1], A[2, 2], A[2, 3]
    a41, a42, a43, a44 = A[3, 0], A[3, 1], A[3, 2], A[3, 3]

    # bij belongs to b = Ma'*(Mb^-1)'*B*(Mb^-1)*Ma .
    aux = np.linalg.inv(Mb) @ Ma
    b = aux.T @ B @ aux
    b11, b12, b13, b14 = b[0, 0], b[0, 1], b[0, 2], b[0, 3]
    b21, b22, b23, b24 = b[1, 0], b[1, 1], b[1, 2], b[1, 3]
    b31, b32, b33, b34 = b[2, 0], b[2, 1], b[2, 2], b[2, 3]
    b41, b42, b43, b44 = b[3, 0], b[3, 1], b[3, 2], b[3, 3]

    # Coefficients of the Characteristic Polynomial.
    T4 = (-a11 * a22 * a33)
    T3 = (a11 * a22 * b33 + a11 * a33 * b22 + a22 * a33 * b11 - a11 * a22 * a33 * b44)
    T2 = (a11 * b23 * b32 - a11 * b22 * b33 - a22 * b11 * b33 + a22 * b13 * b31 -
          a33 * b11 * b22 + a33 * b12 * b21 + a11 * a22 * b33 * b44 - a11 * a22 * b34 * b43 +
          a11 * a33 * b22 * b44 - a11 * a33 * b24 * b42 + a22 * a33 * b11 * b44 -
          a22 * a33 * b14 * b41)
    T1 = (b11 * b22 * b33 - b11 * b23 * b32 - b12 * b21 * b33 + b12 * b23 * b31 +
          b13 * b21 * b32 - b13 * b22 * b31 - a11 * b22 * b33 * b44 + a11 * b22 * b34 * b43 +
          a11 * b23 * b32 * b44 - a11 * b23 * b34 * b42 - a11 * b24 * b32 * b43 +
          a11 * b24 * b33 * b42 - a22 * b11 * b33 * b44 + a22 * b11 * b34 * b43 +
          a22 * b13 * b31 * b44 - a22 * b13 * b34 * b41 - a22 * b14 * b31 * b43 +
          a22 * b14 * b33 * b41 - a33 * b11 * b22 * b44 + a33 * b11 * b24 * b42 +
          a33 * b12 * b21 * b44 - a33 * b12 * b24 * b41 - a33 * b14 * b21 * b42 +
          a33 * b14 * b22 * b41)
    T0 = (b11 * b22 * b33 * b44 - b11 * b22 * b34 * b43 - b11 * b23 * b32 * b44 +
          b11 * b23 * b34 * b42 + b11 * b24 * b32 * b43 - b11 * b24 * b33 * b42 -
          b12 * b21 * b33 * b44 + b12 * b21 * b34 * b43 + b12 * b23 * b31 * b44 -
          b12 * b23 * b34 * b41 - b12 * b24 * b31 * b43 + b12 * b24 * b33 * b41 +
          b13 * b21 * b32 * b44 - b13 * b21 * b34 * b42 - b13 * b22 * b31 * b44 +
          b13 * b22 * b34 * b41 + b13 * b24 * b31 * b42 - b13 * b24 * b32 * b41 -
          b14 * b21 * b32 * b43 + b14 * b21 * b33 * b42 + b14 * b22 * b31 * b43 -
          b14 * b22 * b33 * b41 - b14 * b23 * b31 * b42 + b14 * b23 * b32 * b41)

    # Roots of the characteristic_polynomial (lambda0, ... , lambda4).
    characteristic_polynomial = np.array([T4, T3, T2, T1, T0])
    r = np.roots(characteristic_polynomial)

    # Correct numerical error of real valued polynomial zeros that are
    # accompanied by complex numbers.
    for k in range(4):
        if np.imag(r[k]) <= 10 ** -3:
            r[k] = np.real(r[k])

    # Algebraic condition for contact detection between ellipsoids.
    # Find complex roots.
    complex_roots = np.iscomplex(r)

    # Find the (real) negative roots.
    negative_roots_ids = np.where((~complex_roots) & (r < 0))[0]

    # Contact detection status.
    if len(negative_roots_ids) == 2:
        if r[negative_roots_ids[0]] != r[negative_roots_ids[1]]:
            print('Separation Condition: quadric surfaces are separated.')
            status = 'y'
            return status
        elif np.abs(r[negative_roots_ids[0]] - r[negative_roots_ids[1]]) <= 10 ** -3:
            print('Separation Condition: quadric surfaces share a single contact point.')
            status = 'n'
            return status
    else:
        print('Separation Condition: quadric surfaces are not separated (overlapping).')
        status = 'n'
        return status


if __name__ == '__main__':
    gmsh.initialize()
    d1 = Dot(0, 0, 0)
    d2 = Dot(7, 0, 0)
    angle = Dot(np.pi/4, np.pi/8, 0)
    E1 = Ellipsoid(d1, 1, 4, 4, 4, d1, 1)
    # 2 1 1
    E1.create_mesh()
    E2 = Ellipsoid(d2, 1, 1, 1, 2, angle, 2)
    # 3 2 4
    E2.create_mesh()

    ainv = np.linalg.inv(E1.matrix)
    b = E2.matrix
    m = ainv @ b
    eigvals, eigvecs = np.linalg.eig(m)
    eigvecs_normalized = []
    for v in eigvecs:
        new_v = np.squeeze(np.asarray(v))
        normalized_v = []

        ok = new_v[3] != 0

        if ok:
            normalized_v = [x / new_v[3] for x in new_v]

        ok *= all([x.imag == 0 for x in normalized_v])
        if ok:
            eigvecs_normalized.append(normalized_v)


    print(eigvecs_normalized)
    dots_to_check = []
    for v in eigvecs_normalized:
        coords = [x.imag for x in v]
        d = Dot(coords=coords)
        dots_to_check.append(d)

    for d in dots_to_check:
        print(d_in_ell(d, E1))
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


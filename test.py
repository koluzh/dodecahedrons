import gmsh
import numpy as np
a, b, c = 2, 3, 4
A = b * b * c * c
B = a * a * c * c
C = a * a * b * b
K = - a * a * b * b * c * c
S = np.matrix([[A, 0, 0, 0], [0, B, 0, 0], [0, 0, C, 0], [0, 0, 0, K]])
x, y, z = 0, 0, 0
T = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [x, y, z, 1]])
alpha = 0
beta = np.pi/2
gamma = 0
R_x = np.matrix([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
R_y = np.matrix([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
R_z = np.matrix([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])
R = np.dot(R_x, R_y)
R = np.dot(R, R_z)
temp1 = np.array(R[0, :])
temp1 = temp1.ravel()
temp1 = np.hstack([temp1, [0]])
temp2 = np.array(R[1, :])
temp2 = temp2.ravel()
temp2 = np.hstack([temp2, [0]])
temp3 = np.array(R[2, :])
temp3 = temp3.ravel()
temp3 = np.hstack([temp3, [0]])
R = np.matrix([temp1, temp2, temp3, [0, 0, 0, 1]])
S = np.dot(S, R)
X = np.matrix([0, 0, 0, 1])
X_t = X.transpose()
inside = np.dot(X, S)
inside = np.dot(inside, X_t)

import geom

origin = geom.Dot(0, 0, 0)
angle = geom.Dot(0, 0, 0)
lol_center = geom.Dot(7, 0, 0)

kek = geom.Ellipsoid(origin, 2, 1, 1)
print("matrix:\n", kek.matrix, '\n')
lol = geom.Ellipsoid(lol_center, 3, 2, 4)
print("matrix:\n", lol.matrix, '\n')

ans1, ans2 = geom.ell_intersection(kek, lol)

print(ans1)
print('\n')
print(ans2)

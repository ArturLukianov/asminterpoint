import numpy
import numpy as np
from numpy.linalg import norm

x = np.array([1 / 2, 7 / 2, 1, 2], float)  # first trial solution
alpha = 0.5  # alpha - just constant
A = np.array(
    [[2, 4, 1, 0], [1, 3, 0, -1]], float
)  # coefficients of subject function (constraints)
c = np.array([-1, -1, 0, 0], float)  # coefficients of objective function


i = 1

while True:
    v = x
    D = np.diag(x)  # taking diagonal items

    AA = np.dot(A, D)
    cc = np.dot(D, c)

    I = np.eye(4)

    F = np.dot(AA, np.transpose(AA))
    FI = np.linalg.inv(F)
    H = np.dot(np.transpose(AA), FI)

    P = np.subtract(I, np.dot(H, AA))

    cp = np.dot(P, cc)

    nu = np.absolute(np.min(cp))
    y = np.add(np.ones(4, float), (alpha / nu) * cp)  # x with hat

    yy = np.dot(D, y)  # trial solution for next iteration

    x = yy

    if i == 1 or i == 2 or i == 3 or i == 4:
        print("In iteration  ", i, " we have x = ", x, "\n")
        i = i + 1

    if norm(np.subtract(yy, v), ord=2) < 0.00001:
        break


print("In the last iteration  ", i, "  we have x=  \n", x)

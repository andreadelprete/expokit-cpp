import numpy as np
from balanceMethods import new_balance
from numpy.linalg import norm

np.set_printoptions(linewidth=250, suppress=True)

A = np.array([[-3.80466072e-01,  1.19867645e+06,  8.63170116e-02,
               -2.60920215e+05],
              [1.68822169e+00,  1.56173278e+05,  4.61874965e+05,
               -9.96831380e+05],
              [1.20856994e+00, -8.61700547e-01, -6.34023305e-01,
               -9.35935672e-01],
              [-2.87160379e+05,  4.70987189e+05, -7.37804710e-01,
               6.75828111e-01]])
print(norm(A, 1, axis=0))
print("----")

B, D, Dinv, it = new_balance(A)

print(it)
print(norm(B, 1))
print(B)
print(D)

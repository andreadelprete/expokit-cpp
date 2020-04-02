from random import randrange
import numpy as np
from testingStuff import testStuff
from usefulStuff import computeGains
from balanceMethods import slow_balance
from scipy.linalg import matrix_balance
from numpy.linalg import norm


testStuff()


def testMatrixGammaNegative(A):
    B_scipy, D1 = matrix_balance(A, permute=False)
    B_new, D2, Dinv, it = slow_balance(A)

    norm_scipy = norm(B_scipy, 1)
    norm_new = norm(B_new, 1)

    return computeGains(norm_scipy, norm_new)


N = 4

while True:
    A = np.random.randn(N, N)
    Er = randrange(0, int((N**2)/2))  # How many elemnt we wanna modify
    for r in range(0, Er):
        # Random position
        c = randrange(0, N)
        r = randrange(0, N)

        mul = float(randrange(10e2, 10e5))  # By how much

        A[r, c] *= mul  # Do it

    gamma, squarings_gain = testMatrixGammaNegative(A)

    # Loop until we find the case we need
    if gamma < 0:
        break


print(testMatrixGammaNegative(A))

A = np.array([[3.905994339605555433e+04, 1.423115386082963285e+00, -7.767689234713232145e+05, 2.116088780703028606e+05],
              [1.630857904774933864e-01, 1.485414816171685803e+00, 1.708381725952879515e+00, -8.168065531205611629e-01],
              [-3.732235755447970565e-02, -2.641333416301955032e+10, 2.940791468506479545e-01, 4.523585287019480927e-01],
              [9.687367685924530961e+05, -6.031993146103931114e-01, -2.769428994406575306e-01, 1.509967925770319998e+05]])


print(testMatrixGammaNegative(A))

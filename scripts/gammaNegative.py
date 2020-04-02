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
        pass


testMatrixGammaNegative(A)

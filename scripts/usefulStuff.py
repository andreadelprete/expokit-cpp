import numpy as np
from numpy.linalg import norm
from scipy.linalg import matrix_balance
from balanceMethods import new_balance, slow_balance2  # , slow_balance
from random import randrange


maxnorm = 5.371920351148152


def compute_squarings(l1norm):
    return int(np.ceil(np.log2(l1norm / maxnorm)))


def printStats(M, text):
    l1norm = norm(M, 1)
    s = compute_squarings(l1norm)
#    print(M)
    # print(c / 10e5)
    print('%20s   %10d   %10.2f' % (text, s, np.log2(l1norm)))


def test_matrix(A):
    B_scipy, D1 = matrix_balance(A, permute=False)
    B_new, D2, Dinv, it = slow_balance2(A)
    # B_new, D2, Dinv, it = new_balance(A)

    return computeGains(norm(B_scipy, 1), norm(B_new, 1))


def computeGains(norm1, norm2):
    gamma = (norm1 - norm2) / norm2
    squarings_gain = compute_squarings(norm1) - compute_squarings(norm2)

    return gamma, squarings_gain


def generateStiffMatrix(N):
    A = np.random.randn(N, N)
    Er = randrange(0, int((N**2)/2))  # How many elemnt we wanna modify
    for r in range(0, Er):
        # Random position
        c = randrange(0, N)
        r = randrange(0, N)

        mul = float(randrange(10e2, 10e5))  # By how much

        A[r, c] *= mul  # Do it

    return A

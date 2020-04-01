import numpy as np
from numpy.linalg import norm
from scipy.linalg import matrix_balance
from balanceMethods import new_balance


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
    B_new, D2, Dinv, it = new_balance(A)

    norm_new = norm(B_new, 1)
    norm_scipy = norm(B_scipy, 1)
    gamma = (norm_scipy - norm_new) / norm_new
    squarings_gain = compute_squarings(norm_scipy) - compute_squarings(norm_new)

    return gamma, squarings_gain

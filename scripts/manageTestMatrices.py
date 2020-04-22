import scipy.io
import os
import numpy as np
from random import randrange
from numpy.linalg import norm


# Matrices are returned into an array and are accessed as such:
# array[k]['A']
# also the name of the matrix is present as array[k]['name']

# Disregard previous comment, it return a simple list of matrices
def importMatrices(path):
    matrices = []
    for filename in os.listdir(os.getcwd() + '/' + path):
        m = scipy.io.loadmat(path + '/' + filename)['A'].real.astype('float64')
        # There is the possibility to trace the MatLab name
        # m['name'] = os.path.splitext(filename)[0]  # Strip extension
        matrices.append(m)

    return matrices


# a = importMatrices('scripts/testMatrices')
# print(a[0])


def james2014Generator(N):
    A = np.random.randn(N, N)
    D = np.zeros(N)
    Dinv = np.zeros(N)

    for n in range(N):
        e = randrange(0, 20)
        D[n] = 2**e
        Dinv[n] = 2**-e

    return norm(A, 1), np.diagflat(Dinv) @ A @ np.diagflat(D)


# S = james2014Generator(8)
# print(S)

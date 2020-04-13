from random import randrange
import numpy as np
from testingStuff import testStuff
from usefulStuff import computeGains
from balanceMethods import slow_balance, new_balance
from scipy.linalg import matrix_balance
from numpy.linalg import norm
import multiprocessing


testStuff()


def testMatrixGammaNegative(A):
    B_scipy, D1 = matrix_balance(A, permute=False)
    B_new, D2, Dinv, it = slow_balance(A)

    norm_scipy = norm(B_scipy, 1)
    norm_new = norm(B_new, 1)

    return computeGains(norm_scipy, norm_new)


def testCombinedAlg(A):
    oNorm = norm(A, 1)

    BScipy, D1 = matrix_balance(A, permute=False)
    normScipy = norm(BScipy, 1)

    # Because if the original balance does worst, also the new algorithm is spoiled
    if oNorm < normScipy:
        BNew, D2, D2inv, it = new_balance(A)
    else:
        BNew, D2, D2inv, it = new_balance(BScipy)

    normNew = norm(BNew, 1)

    res = computeGains(normScipy, normNew)
    return (oNorm < normNew), res[0], res[1]


''' Test matrices that lead to local minima solutions
# Very similar approach in both methods
A = np.array([[3.905994339605555433e+04, 1.423115386082963285e+00, -7.767689234713232145e+05, 2.116088780703028606e+05],
              [1.630857904774933864e-01, 1.485414816171685803e+00, 1.708381725952879515e+00, -8.168065531205611629e-01],
              [-3.732235755447970565e-02, -2.641333416301955032e+10, 2.940791468506479545e-01, 4.523585287019480927e-01],
              [9.687367685924530961e+05, -6.031993146103931114e-01, -2.769428994406575306e-01, 1.509967925770319998e+05]])
# print(testMatrixGammaNegative(A))
print(testCombinedAlg(A))

# Completely different approach
A = np.array([[1.194404702374889160e+00, -1.880292268576423521e+04, 5.036505233183925512e-01, -5.910019416967572021e+11],
              [2.303119118211722771e+00, 7.205807592369093406e-01, 5.909011499385113941e-01, 1.323280094717320110e-01],
              [-9.685716355729805382e-01, -1.889628260705311550e+05, 7.994406038426891126e-01, 7.033830337305769786e-01],
              [-1.579883917107374469e-01, 4.153532186150737382e-01, -9.987296987985255781e-01, -7.156691908136625369e-01]])
print(testMatrixGammaNegative(A))
'''

# N = 4


def searchForGammaNegative(N):
    maxIterations = 10e3
    i = 0
    while True:
        A = np.random.randn(N, N)
        Er = randrange(0, int((N**2)/2))  # How many elemnt we wanna modify
        for r in range(0, Er):
            # Random position
            c = randrange(0, N)
            r = randrange(0, N)

            mul = float(randrange(10e2, 10e5))  # By how much

            A[r, c] *= mul  # Do it

        # gamma, squarings_gain = testMatrixGammaNegative(A)
        worse, gamma, squarings_gain = testCombinedAlg(A)

        # Loop until we find the case we need
        if gamma < 0 or worse:
            print("Found one at iteration %i for N = %i" % (i, N), worse, testMatrixGammaNegative(A))
            break

        if i > maxIterations:
            print("Nothing found for N = %i in %i iterations" % (N, maxIterations))
            break

        i = i + 1


# searchForGammaNegative(4)

matrixDimensions = [4, 6, 8, 10, 12, 14]
p = multiprocessing.Pool()
p.map(searchForGammaNegative, matrixDimensions)

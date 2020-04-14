import numpy as np
from testingStuff import testStuff
from usefulStuff import generateStiffMatrix
from balanceMethods import slow_balance, new_balance
from scipy.linalg import matrix_balance
from numpy.linalg import norm

testStuff()


def testBiggerNorm(A):
    oNorm = norm(A, 1)

    BScipy, D1 = matrix_balance(A, permute=False)
    normScipy = norm(BScipy, 1)

    BComb, D2, D2inv, it2 = new_balance(BScipy)
    normComb = norm(BComb, 1)

    BNew, D3, D3inv, it3 = slow_balance(A)
    normNew = norm(BNew, 1)

    test = np.sum(D3) == A.shape[0]

    return oNorm < normScipy and oNorm < normComb # and not test # uncommment for infinite loop


def searchForWorseCases(N):
    while True:
        A = generateStiffMatrix(N)

        if testBiggerNorm(A):
            # return A
            pass


searchForWorseCases(4)
# np.savetxt('worse3.boh', searchForWorseCases(4), delimiter=', ')

# In general all these matrices have one VERY BIG element
# A = np.array([  # Case were if SciPy worsen norm, also NB does
#     [-1.263022039120193507e-01, -8.606493431462458599e-01, -3.419659671159938075e-01, 2.191191911157561456e-01],
#     [-1.301465092615038177e+00, -1.110889986525864960e+17, -6.948565837226625685e-01, 4.461342254067287016e-01],
#     [3.122262686809100796e-01, -3.565990987169408988e-01, 4.645004642195044435e+05, 7.252286542047614930e+05],
#     [-8.697324684124282390e-01, -1.379602591329433192e+00, 3.268784190306180948e-01, 2.558731869222150568e+00]])
# testBiggerNorm(A)

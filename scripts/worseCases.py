import numpy as np
from testingStuff import testStuff
from usefulStuff import generateStiffMatrix
from balanceMethods import new_balance
from scipy.linalg import matrix_balance
from numpy.linalg import norm

np.set_printoptions(precision=2, linewidth=200, suppress=True)
testStuff()


def testBiggerNorm(A):
    oNorm = norm(A, 1)
    print('A\n', A)
    print('norm(A)=', oNorm)

    BNB, D0, D0inv, it0 = new_balance(A)
    normNB = norm(BNB, 1)
    print('\nBNV\n', BNB)
    print('norm(BNV)=', normNB)

    BScipy, D1 = matrix_balance(A, permute=False)
    normScipy = norm(BScipy, 1)
    print('\nBscipy\n', BScipy)
    print('D', np.diag(D1))
    print('norm(Bscipy)=', normScipy)

    BComb, D2, D2inv, it2 = new_balance(BScipy)
    normComb = norm(BComb, 1)
    print('\nBComb\n', BComb)
    print('norm(BComb)=', normComb)

    if(oNorm < normComb):
        gamma = 1.0-normComb/oNorm
        print('gamma Comb =', gamma)
        print('gamma scipy=', 1.0-normScipy/oNorm)
        if(gamma < -0.01):
            return True
    return False


def searchForWorseCases(N):
    while True:
        A = generateStiffMatrix(N)

        if testBiggerNorm(A):
            return A
#            pass


# A = searchForWorseCases(4)
# np.savetxt('worse3.boh', A, delimiter=', ')

# In general all these matrices have one VERY BIG element
A = np.array([  # Case were if SciPy worsen norm, also NB does
    [1.733247591628952478e-01, 6.501614395203411667e-02, 6.181838654801411481e-01, 5.352890587219977236e-01],
    [1.432357938829300403e+00, -4.087364413937133456e-01, -7.975387086316637619e-01, -1.301799804473845690e+00],
    [1.167866623091465383e+00, -1.779488489624770908e+00, 1.732520589971199509e+00, 1.127764779665160422e+00],
    [-8.274426178613810690e-01, -2.877275828268960889e-01, 1.286993711720545719e-01, -6.307435176166014124e-01]])
testBiggerNorm(A)

from random import randrange
import numpy as np
from numpy.linalg import norm
from scipy.linalg import matrix_balance

from testingStuff import testStuff
from balanceMethods import balance_rodney, new_balance, simple_balance


np.set_printoptions(precision=2, linewidth=250, suppress=True)

N = 10
maxnorm = 5.371920351148152

testStuff()


def printStats(M, text):
    l1norm = norm(M, 1)
    s = int(np.ceil(np.log2(l1norm / maxnorm)))
    # print(M)
    # print(c / 10e5)
    print(text, s, "\t%.3f" % l1norm)


for k in range(0, 10):
    A = np.random.randn(N, N)
    Er = randrange(0, (N**2)/2)  # How many elemnt we wanna modify
    for r in range(0, Er):
        # Random position
        c = randrange(0, N)
        r = randrange(0, N)

        mul = float(randrange(10e2, 10e5))  # By how much

        A[r, c] *= mul  # Do it

    printStats(A, "Original\t")

    B, D = matrix_balance(A, permute=False)
    printStats(B, "Scypy\t\t")

    B, D, Dinv, it = balance_rodney(A)
    printStats(B, "Rodney\t\t")

    B, D, Dinv, it = new_balance(A)
    printStats(B, "Our New\t\t")

    B, D, Dinv, it = simple_balance(A)
    printStats(B, "OurSimple\t")

    print("#####################################")

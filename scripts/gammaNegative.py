from random import randrange
import numpy as np
from testingStuff import testStuff
from usefulStuff import test_matrix

testStuff()


N = 4

while True:
    A = np.random.randn(N, N)
    Er = randrange(0, int((N**2)/2))  # How many elemnt we wanna modify
    for r in range(0, Er):
        # Random position
        c = randrange(0, N)
        r = randrange(0, N)

        mul = float(randrange(10e2, 10e3))  # By how much

        A[r, c] *= mul  # Do it
    # Loop until to avoid creation of matrices that does not need balancing

    gamma, squarings_gain = test_matrix(A)

    if gamma < 0:
        break


print(A)

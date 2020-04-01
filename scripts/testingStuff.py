from balanceMethods import balance_rodney, new_balance, simple_balance
import numpy as np
from numpy.linalg import norm
import math


def testStuff():
    # Testing balance methods
    A = np.array([[1, 2, 300], [4, 5, 6], [0, 0, 0]], dtype=np.float)

    B, D, Dinv, it = balance_rodney(A)
    assert(math.isclose(norm(A - D @ B @ Dinv), 0.0, rel_tol=1e-5))

    B, D, Dinv, it = new_balance(A)
    assert(math.isclose(norm(A - D @ B @ Dinv), 0.0, rel_tol=1e-5))

    B, D, Dinv, it = simple_balance(A)
    assert(math.isclose(norm(A - D @ B @ Dinv), 0.0, rel_tol=1e-5))

    print("The stuff you are using makes sense, go on\n")


if __name__ == '__main__':
    testStuff()

import numpy as np
from numpy.linalg import norm


def balance_rodney(A):
    n = A.shape[0]
    assert(A.shape[1] == n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation
    Dinv = np.ones(n)  # Cheaper to compute along the way
    rb = 2.0  # Radix base

    converged = False

    it = 0
    while not converged:
        converged = True
        for i in range(n):
            c = norm(B[:, i], 1)
            r = norm(B[i, :], 1)
            s = c+r
            f = 1.0

            while c < r / rb and c != 0:
                c *= rb
                r /= rb
                f *= rb

            while c > r * rb and r != 0:
                c /= rb
                r *= rb
                f /= rb

            if(c + r < 0.95 * s):
                converged = False
                D[i] *= f
                Dinv[i] /= f
                B[:, i] *= f
                B[i, :] /= f

        it = it + 1

    return B, np.diagflat(D), np.diagflat(Dinv), it


def new_balance(A, max_iter=None):
    n = A.shape[0]
    assert(A.shape[1] == n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation
    Dinv = np.ones(n)  # Cheaper to compute along the way
    rb = 2.0  # Radix base

    # the (i,j) element of cNew is the new 1-norm of column j of the balanced
    # matrix you would get modifying element i of D
    cNew = np.zeros((n, n))

    # each element of v contains the new 1-norm of balanced matrix you would
    # get if you modified the corresponding element of D
    v = np.zeros(n)

    converged = False
    # compute the 1-norm of each column
    c = norm(B, 1, axis=0)
    if(max_iter is None):
        max_iter = n*n

    it = 0
    for it in range(max_iter):
        # find the column with largest norm
        jMax = np.argmax(c)

#        print("\nIter %d |B|=%.1f"%(it, c[jMax]))
#        print("c =", c.T)

        # compute the hypothetical new 1-norms of B
        vMin = 10*c[jMax]
        iMin = 0
        for k in range(n):
            if(k == jMax):
                cNew[jMax, :] = c + np.abs(B[jMax, :])
                cNew[jMax, jMax] = c[jMax] / rb
                v[jMax] = np.max(cNew[jMax, :])
            else:
                cNew[k, :] = c
                cNew[k, k] *= rb
                cNew[k, jMax] -= B[k, jMax] / rb
                v[k] = np.max(cNew[k, :])

            if(v[k] < vMin):
                vMin = v[k]
                iMin = k

#        print("v =", v.T)
#        print("iMin=%d, vMin=%.1f"%(iMin, vMin))

        # check convergence
        if(vMin >= c[jMax]):
            # print("Balancing converged in %d iterations" % it)
            converged = True
            break

        # pick greediest choice
        if(iMin == jMax):
            #            print("Apply strategy 1")
            D[iMin] /= rb
            Dinv[iMin] *= rb
            B[:, iMin] /= rb
            B[iMin, :] *= rb
        else:
            # print("Apply strategy 2")
            D[iMin] *= rb
            Dinv[iMin] /= rb
            B[:, iMin] *= rb
            B[iMin, :] /= rb

        c = np.copy(cNew[iMin, :])

    if(not converged):
        print("ERROR: balancing algorithm did not converge!")

    return B, np.diagflat(D), np.diagflat(Dinv), it

import numpy as np
from numpy.linalg import norm
import math


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
    jMax = np.argmax(c)
    if(max_iter is None):
        max_iter = 2+n*int(2+np.log2(c[jMax]))

    it = 0
#    print('A\n', A)
    for it in range(max_iter):
        # find the column with largest norm
        #        c = norm(B, 1, axis=0)

        #        print("\nIter %d |B|=%.1f, jMax=%d"%(it, c[jMax], jMax))
        #        print('B[:,jMax] = ', B[:,jMax].T)
        #        print("c =", c.T)

        # compute the hypothetical new 1-norms of B
        vMin = 10*c[jMax]
        iMin = 0
        for k in range(n):
            if(k == jMax):
                # the elements of row jMax are multiplied by rb
                cNew[k, :] = c + (rb-1)*np.abs(B[k, :])
                # all elements of col jMax are divided by rb (except for the one on the diagonal)
                cNew[k, k] = (c[k] + (rb-1)*abs(A[k, k])) / rb
                v[k] = np.max(cNew[k, :])
            else:
                # the elements of row k are divided by rb
                cNew[k, :] = c - (rb-1)*np.abs(B[k, :])/rb
                # the elements of column k are multiplied by rb (except for the one on the diagonal)
                cNew[k, k] = c[k]*rb - (rb-1)*abs(A[k, k])
                v[k] = np.max(cNew[k, :])

            if(v[k] < vMin):
                vMin = v[k]
                iMin = k
#        print('cNew = \n', cNew)
#        print("v =", v.T)
#        print("iMin=%d, vMin=%.1f"%(iMin, vMin))

        # check convergence
        if(vMin >= c[jMax]):
            #            print("Balancing converged in %d iterations" % it)
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
            #            print("Apply strategy 2")
            D[iMin] *= rb
            Dinv[iMin] /= rb
            B[:, iMin] *= rb
            B[iMin, :] /= rb

#        print('cNew[iMin] = ', cNew[iMin,:])
#        print('cReal =      ', norm(B, 1, axis=0))
        c = np.copy(cNew[iMin, :])
        jMax = np.argmax(c)

    if(not converged):
        print("ERROR: balancing algorithm did not converge!")

    return B, np.diagflat(D), np.diagflat(Dinv), it


def simple_balance(A):
    n = A.shape[0]
    assert(A.shape[1] == n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation
    Dinv = np.ones(n)  # Cheaper to compute along the way
    rb = 2.0  # Radix base

    prevNorm = norm(B, 1)

    it = 0
    warning = 0
    while True:
        it = it + 1

        cAll = norm(B, 1, axis=0)  # Norm of all columns
        maxI = np.argmax(cAll)  # Which column is the bigger
        # Applying changes
        D[maxI] /= rb
        Dinv[maxI] *= rb
        B[:, maxI] /= rb
        B[maxI, :] *= rb

        # Check if the norm of the matrix is better
        nowNorm = norm(B, 1)

        if math.isclose(nowNorm, prevNorm, rel_tol=0.4):
            warning = warning + 1
        else:
            warning = 0

        prevNorm = nowNorm

        if warning > 5:
            break  # And get out of here

    return B, np.diagflat(D), np.diagflat(Dinv), it


# Simplest possibly slowest way of writing new balance logic
def slow_balance(A):
    n = A.shape[0]
    assert(A.shape[1] == n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation
    Dinv = np.ones(n)  # Cheaper to compute along the way
    rb = 2.0  # Radix base

    v1 = np.zeros(n)
    v2 = np.zeros(n)

    it = 0
    while True:
        origNorm = norm(B, 1)

        # Try every possible change to the transformation matrix
        for i in range(n):
            DTemp = np.eye(n)
            DinvTemp = np.eye(n)
            DTemp[i, i] = rb
            DinvTemp[i, i] = 1 / rb

            # Strategy 1
            v1[i] = norm(DTemp @ B @ DinvTemp, 1)

            # Strategy 2
            v2[i] = norm(DinvTemp @ B @ DTemp, 1)

        # Choose the best
        strategy = bool(np.argmin(np.array([v1.min(), v2.min()])))

        if origNorm > min([v1.min(), v2.min()]):  # Only if it leads to an improvement
            if not strategy:  # Thus Strategy 1
                k = np.argmin(v1)
                D[k] /= rb
                Dinv[k] *= rb
                B[:, k] /= rb
                B[k, :] *= rb
            else:  # Strategy 2
                k = np.argmin(v2)
                D[k] *= rb
                Dinv[k] /= rb
                B[:, k] *= rb
                B[k, :] /= rb
        else:
            break

        it = it + 1

    return B, np.diagflat(D), np.diagflat(Dinv), it

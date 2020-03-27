import numpy as np
from numpy.linalg import norm


def balance_parlett(A, max_iter=None):
    n = A.shape[0]
    assert(A.shape[1] == n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation
    converged = False
    if(max_iter is None):
        max_iter = n*n

    for it in range(max_iter):
        converged = True
        for i in range(n):
            c = norm(B[:, i], 1)
            r = norm(B[i, :], 1)
            s = c+r
            f = 1

            while c < 0.5*r:
                c *= 2.0
                r *= 0.5
                f *= 2.0
                if(f > 2**20):
                    break
            while c > 2.0*r:
                c *= 0.5
                r *= 2.0
                f *= 0.5
                if(f < 2**(-20)):
                    break

            if(c+r < 0.95*s):
                if(D[i]*f < 2**(-20)):
                    f = 2**(-20) / D[i]
                    if(f != 1.0):
                        converged = False
                elif(D[i]*f > 2**(20)):
                    f = 2**(20) / D[i]
                    if(f != 1.0):
                        converged = False
                else:
                    converged = False
                D[i] *= f
                B[:, i] *= f
                B[i, :] /= f

        if(converged):
            print("Balancing converged in %d iterations" % it)
            break

    if(not converged):
        print("ERROR: balancing algorithm did not converge!")

    # check that B = Dinv * A * D
#    Dm = np.diagflat(D)
#    Dinv = np.linalg.inv(Dm)
#    B2 = Dinv * A * Dm
#    print('A\n', A)
    print('log2(D)\n', np.log2(D).T)
#    print('B\n', B)
#    print('Dinv*A*D\n', B2)
    return B, np.diagflat(D)


def new_balance(A, max_iter=None):
    n = A.shape[0]
    assert(A.shape[1] == n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation

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
                cNew[jMax, :] = c+np.abs(B[jMax, :])
                cNew[jMax, jMax] = 0.5*c[jMax]
                v[jMax] = np.max(cNew[jMax, :])
            else:
                cNew[k, :] = c
                cNew[k, k] *= 2.0
                cNew[k, jMax] -= 0.5*B[k, jMax]
                v[k] = np.max(cNew[k, :])

            if(v[k] < vMin):
                vMin = v[k]
                iMin = k

#        print("v =", v.T)
#        print("iMin=%d, vMin=%.1f"%(iMin, vMin))

        # check convergence
        if(vMin >= c[jMax]):
            print("Balancing converged in %d iterations" % it)
            converged = True
            break

        # pick greediest choice
        if(iMin == jMax):
            #            print("Apply strategy 1")
            D[iMin] *= 0.5
            B[:, iMin] *= 0.5
            B[iMin, :] *= 2
        else:
            print("Apply strategy 2")
            D[iMin] *= 2
            B[:, iMin] *= 2
            B[iMin, :] *= 0.5

        c = np.copy(cNew[iMin, :])

    if(not converged):
        print("ERROR: balancing algorithm did not converge!")

    # check that B = Dinv * A * D
#    Dm = np.diagflat(D)
#    Dinv = np.linalg.inv(Dm)
#    B2 = Dinv * A * Dm
#    print('A\n', A)
    print('log2(D)\n', np.log2(D).T)
#    print('B\n', B)
#    print('Dinv*A*D\n', B2)
    return B, np.diagflat(D)

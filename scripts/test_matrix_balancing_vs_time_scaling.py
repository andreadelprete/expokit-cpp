'''
Consider a Linear Dynamical System:
  dx(t)/dt = A * x(t) + a
Given x(0) and T I want to compute:
- x(T)
- The integral of x(t) for t in [0, T]
- The double integral of x(t) for t in [0, T]
'''
from __future__ import print_function
from scipy.sparse.linalg.matfuncs import _ExpmPadeHelper, _ell, _solve_P_Q, _onenorm
import numpy as np
from numpy.linalg import norm
from numpy import matlib
from numpy.linalg import solve
from scipy.linalg import matrix_balance
# rom scipy.linalg import expm
from numpy.linalg import eigvals
# import matplotlib.pyplot as plt
#from utils_LDS_integral import compute_x_T
from math import sqrt

np.set_printoptions(precision=2, linewidth=200, suppress=True)

def expm(A, use_exact_onenorm="auto", verbose=False):
    # Core of expm, separated to allow testing exact and approximate
    # algorithms.
    # Hardcode a matrix order threshold for exact vs. estimated one-norms.
    use_exact_onenorm = A.shape[0] < 200
    h = _ExpmPadeHelper(A, use_exact_onenorm=use_exact_onenorm)
    # Use Pade order 13.
    l1norm = _onenorm(A)
    theta_13 = 4.25
    # Choose smallest s>=0 such that 2**(-s) l1norm <= theta_13
    s = max(int(np.ceil(np.log2(l1norm / theta_13))), 0)
#    s = s + _ell(2**-s * h.A, 13)
    print("Number of squarings", s, ". L-1 norm %.3f"%l1norm)
    U, V = h.pade13_scaled(s)
    X = _solve_P_Q(U, V)
    # X = r_13(A)^(2^s) by repeated squaring.
    for i in range(s):
        X = X.dot(X)
    return X

def balance_parlett(A, max_iter=None):
    n = A.shape[0]
    assert(A.shape[1]==n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation    
    converged = False
    if(max_iter is None):
        max_iter = n*n
        
    for it in range(max_iter):        
        converged = True
        for i in range(n):
            c = norm(B[:,i], 1)
            r = norm(B[i,:], 1)
            s = c+r
            f = 1

            while c < 0.5*r:
                c *= 2.0
                r *= 0.5
                f *= 2.0
                if(f>2**20):
                    break
            while c > 2.0*r:
                c *= 0.5
                r *= 2.0
                f *= 0.5
                if(f<2**(-20)):
                    break
            
            if(c+r < 0.95*s):
                if(D[i]*f<2**(-20)):
                    f = 2**(-20) / D[i]
                    if(f!=1.0):
                        converged = False
                elif(D[i]*f>2**(20)):
                    f = 2**(20) / D[i]
                    if(f!=1.0):
                        converged = False
                else:
                    converged = False
                D[i] *= f
                B[:,i] *= f
                B[i,:] /= f
        
        if(converged):
            print("Balancing converged in %d iterations"%it)    
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
    assert(A.shape[1]==n)
    B = np.copy(A)  # balanced matrix
    D = np.ones(n)  # diagonal elements of similarity transformation

    # the (i,j) element of cNew is the new 1-norm of column j of the balanced 
    # matrix you would get modifying element i of D
    cNew = np.zeros((n,n))    

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
            if(k==jMax):
                cNew[jMax, :] = c+np.abs(B[jMax,:])
                cNew[jMax, jMax] = 0.5*c[jMax]
                v[jMax] = np.max(cNew[jMax,:])
            else:
                cNew[k,:] = c
                cNew[k,k] *= 2.0
                cNew[k,jMax] -= 0.5*B[k,jMax]
                v[k] = np.max(cNew[k,:])
            
            if(v[k] < vMin):
                vMin = v[k]
                iMin = k
        
#        print("v =", v.T)
#        print("iMin=%d, vMin=%.1f"%(iMin, vMin))
        
        # check convergence
        if(vMin >= c[jMax]):
            print("Balancing converged in %d iterations"%it)
            converged = True
            break
        
        # pick greediest choice
        if(iMin==jMax):
#            print("Apply strategy 1")
            D[iMin] *= 0.5
            B[:,iMin] *= 0.5
            B[iMin,:] *= 2
        else:
            print("Apply strategy 2")
            D[iMin] *= 2
            B[:,iMin] *= 2
            B[iMin,:] *= 0.5
            
        c = np.copy(cNew[iMin,:])
            
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
            

def print_error(x_exact, x_approx):
    print("Approximation error: ", np.max(np.abs(x_exact-x_approx).A1 / np.abs(x_exact).A1))
    
def compute_x_T(A, a, x0, T, balance=None):
    n = A.shape[0]
    C = matlib.zeros((n+1, n+1))
    C[0:n,     0:n] = A
    C[0:n,     n] = a
    z0 = matlib.zeros((n+1, 1))
    z0[:n, 0] = x0
    z0[-1, 0] = 1.0
    if balance:
        if(balance=='lapack'):
#            C_bal, D = matrix_balance(T*C, permute=False)
            C_bal, D = balance_parlett(T*C)
        else:
            (C_bal, D) = new_balance(T*C)
        D = np.asmatrix(D)
        C_bal = np.asmatrix(C_bal)
        Dinv = np.asmatrix(np.linalg.inv(D))
            
        e_TC_bal = expm(C_bal, verbose=True)
#        print("A\n", T*C)
#        print("Ab\n", C_bal)
#        print("log(D)\n", np.log(np.diag(D)))
        e_TC = D*e_TC_bal*Dinv
    else:
        e_TC = expm(T*C, verbose=True)
    
    z = e_TC*z0
    x_T = z[:n, 0]
    return x_T
    
class LD2S_Util:
    
    def __init__(self):
        pass
    
    def setScaling(self, K):
        l1norm = _onenorm(K)
        self.Ts = 2**int(np.log2(sqrt(1. / l1norm)))

        # Structure to scale only the latter half (speed) of xInit
        n = K.shape[0]
        self.mulXInit = matlib.ones((n*2,1))
        self.mulXInit[n:] *= self.Ts
        
    def composeA(self, K, B):
        n = self.n
        self.As = matlib.zeros((2*n, 2*n))
        self.As[:n, n:] = matlib.identity(n)
        self.As[n:, n:] = -B * self.Ts
        self.As[n:, :n] = -K * (self.Ts**2)
    
    def compute_x_T(self, K, B, b, x0, T):
        self.n = K.shape[0]
        self.setScaling(K)
        self.composeA(K, B)
        x0_s = np.multiply(x0, self.mulXInit)
        out = compute_x_T(self.As, b*self.Ts, x0_s, T/self.Ts)
        out[self.n:] /= self.Ts
        return out


if __name__ == '__main__':
    import time
    N_TESTS = 1
    T = 0.01
    n = 4*3*2
    n2 = int(n/2)
    stiffness = 1e5
    damping = 1e2
    x0 = matlib.rand((n, 1))
    a = matlib.rand((n, 1))
    a[:n2] = 0.0
#    a[:] = 0.0      # TEMP
    U = matlib.rand((n2, n2))
    Upsilon = U*U.T
    K = matlib.eye(n2)*stiffness
    B = matlib.eye(n2)*damping
    A = matlib.block([[matlib.zeros((n2, n2)), matlib.eye(n2)],
                      [-Upsilon*K,      -Upsilon*B]])
    # A  = matlib.rand((n, n))
    util = LD2S_Util()

    print("x(0) is:", x0.T)
    print("a is:   ", a.T)
    print("State size n:", n)
    print("Eigenvalues of A:", np.sort_complex(eigvals(A)).T)
    print("")

#    start_time = time.time()
#    e_TA = expm(T*A)
#    time_exp = time.time()-start_time
#    print("Time to compute matrix exponential", 1e3*time_exp)
#
#    start_time = time.time()
#    A_inv_a = solve(A, a)
#    time_solve = time.time()-start_time
#    print("Time to solve linear system", 1e3*time_solve)
    print("")
    
#    (A_bal, D) = new_balance(A)
    print("Compute expm with new balancing")
    x_T_bal_new = compute_x_T(A, a, x0, T, balance='new')
    
    print("\nCompute expm with standard balancing")
    x_T_bal = compute_x_T(A, a, x0, T, balance='lapack')

    Kbar = Upsilon*K
    Bbar = Upsilon*B
    print("\nCompute expm with time scaling")
    start_time = time.time()
    x_T_ts = util.compute_x_T(Kbar, Bbar, a, x0, T)
    time_approx = time.time()-start_time

    print("\nCompute expm without any balancing")
    start_time = time.time()
    x_T = compute_x_T(A, a, x0, T)
    time_exact = time.time()-start_time
    
    print("\nTime-scaled x(T) computed in             ", 1e3*time_approx)
    print("Standard x(T) computed in                ", 1e3*time_exact)
    print_error(x_T, x_T_ts)
    print("")
    print("x standard  ", x_T.T)
    print("x TimeScaled", x_T_ts.T)
    print("x Balancing ", x_T_bal.T)
    print("x Balancing ", x_T_bal_new.T)

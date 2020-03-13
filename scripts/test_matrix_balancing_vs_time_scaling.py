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
    print("Number of squarings", s, ". L-1 norm", l1norm)
    U, V = h.pade13_scaled(s)
    X = _solve_P_Q(U, V)
    # X = r_13(A)^(2^s) by repeated squaring.
    for i in range(s):
        X = X.dot(X)
    return X


def print_error(x_exact, x_approx):
    print("Approximation error: ", np.max(np.abs(x_exact-x_approx).A1 / np.abs(x_exact).A1))
    
def compute_x_T(A, a, x0, T, balance=False):
    n = A.shape[0]
    C = matlib.zeros((n+1, n+1))
    C[0:n,     0:n] = A
    C[0:n,     n] = a
    z0 = matlib.zeros((n+1, 1))
    z0[:n, 0] = x0
    z0[-1, 0] = 1.0
    if balance:
        C_bal, D = matrix_balance(T*C, permute=False)
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
        self.Ts = sqrt(1. / l1norm);

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
    n = 1*3*2
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

    start_time = time.time()
    e_TA = expm(T*A)
    time_exp = time.time()-start_time
    print("Time to compute matrix exponential", 1e3*time_exp)

    start_time = time.time()
    A_inv_a = solve(A, a)
    time_solve = time.time()-start_time
    print("Time to solve linear system", 1e3*time_solve)
    print("")
    
    x_T_bal = compute_x_T(A, a, x0, T, balance=True)

    Kbar = Upsilon*K
    Bbar = Upsilon*B
    start_time = time.time()
    x_T_ts = util.compute_x_T(Kbar, Bbar, a, x0, T)
    time_approx = time.time()-start_time

    start_time = time.time()
    x_T = compute_x_T(A, a, x0, T)
    time_exact = time.time()-start_time
    print("Time-scaled x(T) computed in             ", 1e3*time_approx)
    print("Standard x(T) computed in                ", 1e3*time_exact)
    print_error(x_T, x_T_ts)
    print("")
    print("x standard  ", x_T.T)
    print("x TimeScaled", x_T_ts.T)
    print("x Balancing ", x_T_bal.T)

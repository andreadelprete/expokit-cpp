import numpy as np
from scipy.sparse.linalg.matfuncs import _ExpmPadeHelper, _solve_P_Q, _onenorm
from math import sqrt
from numpy import matlib
from balanceMethods import balance_rodney, new_balance


class LDSUtil:

    def __init__(self):
        pass

    def setScaling(self, K):
        l1norm = _onenorm(K)
        self.Ts = 2**int(np.log2(sqrt(1. / l1norm)))

        # Structure to scale only the latter half (speed) of xInit
        n = K.shape[0]
        self.mulXInit = matlib.ones((n*2, 1))
        self.mulXInit[n:] *= self.Ts

    def composeA(self, K, B):
        n = self.n
        self.As = matlib.zeros((2*n, 2*n))
        self.As[:n, n:] = matlib.identity(n)
        self.As[n:, n:] = -B * self.Ts
        self.As[n:, :n] = -K * (self.Ts**2)

    def computeXTTimeScaling(self, K, B, b, x0, T):
        self.n = K.shape[0]
        self.setScaling(K)
        self.composeA(K, B)
        x0_s = np.multiply(x0, self.mulXInit)
        out = self.computeXT(self.As, b*self.Ts, x0_s, T/self.Ts)
        out[self.n:] /= self.Ts
        return out

    def computeXT(self, A, a, x0, T, balance=None):
        n = A.shape[0]
        C = matlib.zeros((n+1, n+1))
        C[0:n,     0:n] = A
        C[0:n,     n] = a
        z0 = matlib.zeros((n+1, 1))
        z0[:n, 0] = x0
        z0[-1, 0] = 1.0
        if balance:
            if(balance == 'lapack'):
                #            C_bal, D = matrix_balance(T*C, permute=False)
                C_bal, D = balance_rodney(T*C)
            else:
                (C_bal, D) = new_balance(T*C)
            D = np.asmatrix(D)
            C_bal = np.asmatrix(C_bal)
            Dinv = np.asmatrix(np.linalg.inv(D))

            e_TC_bal = self.expm(C_bal, verbose=True)
    #        print("A\n", T*C)
    #        print("Ab\n", C_bal)
    #        print("log(D)\n", np.log(np.diag(D)))
            e_TC = D*e_TC_bal*Dinv
        else:
            e_TC = self.expm(T*C, verbose=True)

        z = e_TC*z0
        x_T = z[:n, 0]
        return x_T

    def expm(self, A, use_exact_onenorm="auto", verbose=False):
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
        print("Number of squarings", s, ". L-1 norm %.3f" % l1norm)
        U, V = h.pade13_scaled(s)
        X = _solve_P_Q(U, V)
        # X = r_13(A)^(2^s) by repeated squaring.
        for i in range(s):
            X = X.dot(X)
        return X

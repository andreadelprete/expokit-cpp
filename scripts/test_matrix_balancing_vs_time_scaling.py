'''
Consider a Linear Dynamical System:
  dx(t)/dt = A * x(t) + a
Given x(0) and T I want to compute:
- x(T)
- The integral of x(t) for t in [0, T]
- The double integral of x(t) for t in [0, T]
'''
import numpy as np
from numpy import matlib
from numpy.linalg import eigvals
import time

from LDSUtil import LDSUtil
from testingStuff import testStuff


np.set_printoptions(precision=2, linewidth=200, suppress=True)


def print_error(x_exact, x_approx):
    print("Approximation error: ", np.max(np.abs(x_exact-x_approx).A1 / np.abs(x_exact).A1))


#testStuff()  # Avoid regression

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
util = LDSUtil()

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
x_T_bal_new = util.computeXT(A, a, x0, T, balance='new')

print("\nCompute expm with standard balancing")
x_T_bal = util.computeXT(A, a, x0, T, balance='lapack')

Kbar = Upsilon*K
Bbar = Upsilon*B
print("\nCompute expm with time scaling")
start_time = time.time()
x_T_ts = util.computeXTTimeScaling(Kbar, Bbar, a, x0, T)
time_approx = time.time()-start_time

print("\nCompute expm without any balancing")
start_time = time.time()
x_T = util.computeXT(A, a, x0, T)
time_exact = time.time()-start_time

print("\nTime-scaled x(T) computed in             ", 1e3*time_approx)
print("Standard x(T) computed in                ", 1e3*time_exact)
print_error(x_T, x_T_ts)
print("")
print("x standard        ", x_T.T)
print("x TimeScaled      ", x_T_ts.T)
print("x Balancing LAPACK", x_T_bal.T)
print("x Balancing NEW   ", x_T_bal_new.T)

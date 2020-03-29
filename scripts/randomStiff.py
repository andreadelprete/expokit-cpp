from random import randrange
import numpy as np
from numpy.linalg import norm
from scipy.linalg import matrix_balance
import plot_utils as plut
import matplotlib.pyplot as plt
from testingStuff import testStuff
from balanceMethods import balance_rodney, new_balance, simple_balance


np.set_printoptions(precision=2, linewidth=250, suppress=True)

MATRIX_SIZES = [2, 3, 4, 10]          # matrix size
N_TESTS = 1000    # number of tests
maxnorm = 5.371920351148152

testStuff()

def compute_squarings(l1norm):
    return int(np.ceil(np.log2(l1norm / maxnorm)))

def printStats(M, text):
    l1norm = norm(M, 1)
    s = compute_squarings(l1norm)
#    print(M)
    # print(c / 10e5)
    print('%20s   %10d   %10.2f' % (text, s, np.log2(l1norm)))

def test_matrix(A):
#    print('%20s | %10s | %10s' % ('Method', 'Squarings', 'log2(L1-norm)'))
#    printStats(A, "Original")

    B_scipy, D = matrix_balance(A, permute=False)
#    printStats(B_scipy, "Scypy")

#    B_rodney, D, Dinv, it = balance_rodney(A)
#    printStats(B_rodney, "Rodney")

    B_new, D, Dinv, it = new_balance(A)
#    printStats(B_new, "Our New")

#    try:
#        B, D, Dinv, it = simple_balance(A)
#        printStats(B, "OurSimple")
#    except:
#        print('OurSimple FAILED')
    
    norm_new = norm(B_new,1)
    norm_scipy = norm(B_scipy,1)
    gamma = (norm_scipy - norm_new) / norm_new
    squarings_gain = compute_squarings(norm_scipy) - compute_squarings(norm_new)
#    print('Gamma = %.3f'%gamma)
#    print("#############################################################")
    return gamma, squarings_gain

def run_random_tests(matrix_size, n_tests):
    N = matrix_size
    gamma = np.empty(n_tests)
    squarings_gain = np.empty(n_tests)

    for k in range(0, n_tests):
        if(k%100==0):
            print('Test', k)
        A = np.random.randn(N, N)
        Er = randrange(0, int((N**2)/2))  # How many elemnt we wanna modify
        for r in range(0, Er):
            # Random position
            c = randrange(0, N)
            r = randrange(0, N)
    
            mul = float(randrange(10e2, 10e5))  # By how much
    
            A[r, c] *= mul  # Do it
    
        gamma[k], squarings_gain[k] = test_matrix(A)
        
    print("#############################################################")
    print('Matrix size = ', N)
    print('Gamma min-avg-max = %7.3f  %7.3f  %7.3f'%(np.min(gamma), np.mean(gamma), np.max(gamma)))
    print('Percentage of negative gamma = ', 1e2*np.count_nonzero(gamma<0)/n_tests)
    print('Squarings gain min-avg-max = %7.3f  %7.3f  %7.3f'%(np.min(squarings_gain), 
                                                              np.mean(squarings_gain), 
                                                              np.max(squarings_gain)))
    print('Percentage of negative squarings gain = ', 1e2*np.count_nonzero(squarings_gain<0)/n_tests)
    
    # the histogram of the data
    plt.figure()
    n, bins, patches = plt.hist(gamma, 30, facecolor='g', alpha=0.75)
    plt.xlabel('Gamma')
    plt.ylabel('Frequency')
    plt.title('Relative gain in L1 norm for size %d'%N)
    plt.grid(True)
    
    plt.figure()
    n, bins, patches = plt.hist(squarings_gain, 20, facecolor='g', alpha=0.75)
    plt.xlabel('Squarings gain')
    plt.ylabel('Frequency')
    plt.title('Number of squarings gained for size %d'%N)
    plt.grid(True)
    
# TEST CASE STUDY 2b FROM PDF
#p = 10
#A = np.array([[1.5*2**p, 2**p],
#              [    2**p,  0.0]])
#test_matrix(A)


# RANDOM TESTS
np.random.seed(19680801) # Fixing random state for reproducibility
for n in MATRIX_SIZES:
    run_random_tests(n, N_TESTS)
plt.show()

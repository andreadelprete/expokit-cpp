import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from testingStuff import testStuff
from usefulStuff import test_matrix, maxnorm, generateStiffMatrix


np.set_printoptions(precision=2, linewidth=250, suppress=True)

MATRIX_SIZES = [2, 3, 4, 10]          # matrix size
N_TESTS = 1000    # number of tests

testStuff()


def run_random_tests(matrix_size, n_tests):
    N = matrix_size
    gamma = np.empty(n_tests)
    squarings_gain = np.empty(n_tests)

    for k in range(0, n_tests):
        if(k % 100 == 0):
            print('Test', k)

        while True:
            A = generateStiffMatrix(N)

            # Loop until to avoid creation of matrices that does not need balancing
            if norm(A, 1) > maxnorm:
                break

        gamma[k], squarings_gain[k] = test_matrix(A)
        if gamma[k] < 0:
            test_matrix(A)
            pass

    print('Matrix size = ', N)
    print('Gamma min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(gamma), np.mean(gamma), np.max(gamma)))
    print('Percentage of negative gamma = ', 1e2*np.count_nonzero(gamma < 0)/n_tests)
    print('Squarings gain min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(squarings_gain),
                                                                np.mean(squarings_gain),
                                                                np.max(squarings_gain)))
    print('Percentage of negative squarings gain = ', 1e2*np.count_nonzero(squarings_gain < 0)/n_tests)
    print("#############################################################\n")

    # the histogram of the data
    plt.figure()
    n, bins, patches = plt.hist(gamma, 30, facecolor='g', alpha=0.75)
    plt.xlabel('Gamma')
    plt.ylabel('Frequency')
    plt.title('Relative gain in L1 norm for size %d' % N)
    plt.grid(True)

    plt.figure()
    n, bins, patches = plt.hist(squarings_gain, 20, facecolor='g', alpha=0.75)
    plt.xlabel('Squarings gain')
    plt.ylabel('Frequency')
    plt.title('Number of squarings gained for size %d' % N)
    plt.grid(True)


# RANDOM TESTS
np.random.seed(19680801)  # Fixing random state for reproducibility
for n in MATRIX_SIZES:
    run_random_tests(n, N_TESTS)
plt.show()

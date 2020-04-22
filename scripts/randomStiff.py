import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from testingStuff import testStuff
from usefulStuff import test_matrix, maxnorm  # , generateStiffMatrix
from manageTestMatrices import james2014Generator, importMatrices
import multiprocessing
from functools import partial


np.set_printoptions(precision=2, linewidth=250, suppress=True)


testStuff()


def run_random_tests(matrix_size, n_tests):
    N = matrix_size
    gamma = np.empty(n_tests)
    groundTs = np.empty(n_tests)
    squarings_gain = np.empty(n_tests)

    for k in range(0, n_tests):
        if(k % 100 == 0):
            print('Test', k)

        while True:
            gt, A = james2014Generator(N)  # generateStiffMatrix(N)

            # Loop until to avoid creation of matrices that does not need balancing
            if norm(A, 1) > maxnorm:
                break

        gamma[k], squarings_gain[k], groundTs[k] = test_matrix(A, gt)  # Be aware of what method is called here

    printData(gamma, squarings_gain, groundTs)


def analyseTheseMatrices(matrices):
    howMany = len(matrices)
    gamma = np.empty(howMany)
    squarings_gain = np.empty(howMany)

    for k in range(0, howMany):
        gamma[k], squarings_gain[k] = test_matrix(matrices[k])

    printData(gamma, squarings_gain, None, howMany, matrices[0].shape[0])


def printData(gamma, squarings_gain, groundTs, n_tests, N):
    print('Matrix size = ', N)
    print('Gamma min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(gamma), np.mean(gamma), np.max(gamma)))
    
    if groundTs is not None:
        print('groundTs min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(groundTs), np.mean(groundTs), np.max(groundTs)))
        print('Percentage of negative gamma = ', 1e2*np.count_nonzero(gamma < 0)/n_tests)
        print('Percentage of negative groundTs = ', 1e2*np.count_nonzero(groundTs < 0)/n_tests)
    
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

    if (groundTs is not None):
        plt.figure()
        n, bins, patches = plt.hist(groundTs, 30, facecolor='g', alpha=0.75)
        plt.xlabel('Ground Truth')
        plt.ylabel('Frequency')
        plt.title('Difference from original norm %d' % N)
        plt.grid(True)

    plt.figure()
    n, bins, patches = plt.hist(squarings_gain, 20, facecolor='g', alpha=0.75)
    plt.xlabel('Squarings gain')
    plt.ylabel('Frequency')
    plt.title('Number of squarings gained for size %d' % N)
    plt.grid(True)

    plt.show()


whatToDo = False

if whatToDo:
    # Section with tests on generated matrices
    MATRIX_SIZES = [4, 8, 16, 32, 64, 128, 256]  # , 512, 1024]  # matrix size
    N_TESTS = 10  # number of tests

    np.random.seed(19680801)  # Fixing random state for reproducibility

    p = multiprocessing.Pool()
    p.map(partial(run_random_tests, n_tests=N_TESTS), MATRIX_SIZES)

    # for n in MATRIX_SIZES:
    #     run_random_tests(n, N_TESTS)
    # plt.show()
else:
    # Section with tests on Toolbox matrices
    matrices = importMatrices('scripts/testMatrices')

    analyseTheseMatrices(matrices)

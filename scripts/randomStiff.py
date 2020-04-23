import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from testingStuff import testStuff
from usefulStuff import test_matrix, test_matrix_GT, maxnorm  # , generateStiffMatrix
from manageTestMatrices import james2014Generator, importMatrices
import multiprocessing
from functools import partial
from balanceMethods import new_balance, combinedOldNew
import psutil  # Used to get number of cores, physical instead of logical


np.set_printoptions(precision=2, linewidth=250, suppress=True)


testStuff()


def run_random_tests(matrix_size, n_tests=None):
    N = matrix_size
    # Do less test the bigger the matrix is
    if n_tests is None:
        n_tests = 10240 // N
    gamma = np.empty((n_tests, 2))
    groundTs = np.empty((n_tests, 2))
    squarings_gain = np.empty((n_tests, 2))
    squarings_gain_GT = np.empty((n_tests, 2))

    for k in range(0, n_tests):
        if(k % 500 == 0 and k is not 0):
            print('Test {} for dimension {}'.format(k, N))

        while True:
            gt, A = james2014Generator(N)  # generateStiffMatrix(N)

            # Loop until to avoid creation of matrices that does not need balancing
            if norm(A, 1) > maxnorm:
                break

        # Be aware of what method is called here
        gamma[k][0], squarings_gain[k][0], groundTs[k][0], squarings_gain_GT[k][0] = test_matrix_GT(A, gt, new_balance)
        gamma[k][1], squarings_gain[k][1], groundTs[k][1], squarings_gain_GT[k][1] = test_matrix_GT(A, gt, combinedOldNew)

    printData("NB", gamma[:, 0], squarings_gain[:, 0], groundTs[:, 0], squarings_gain_GT[:, 0], n_tests, N)
    printData("NBComb", gamma[:, 1], squarings_gain[:, 1], groundTs[:, 1], squarings_gain_GT[:, 1], n_tests, N)


def analyseTheseMatrices(matrices):
    howMany = len(matrices)
    gamma = np.empty(howMany)
    squarings_gain = np.empty(howMany)

    for k in range(0, howMany):
        gamma[k], squarings_gain[k] = test_matrix(matrices[k])

    printData("Toolbox Matrices", gamma, squarings_gain, None, None, howMany, matrices[0].shape[0])


def printData(title, gamma, squarings_gain, groundTs, squarings_gain_GT, n_tests, N):
    print('Matrix size = ', N)
    print('Gamma min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(gamma), np.mean(gamma), np.max(gamma)))

    if groundTs is not None and squarings_gain_GT is not None:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
        fig.set_size_inches(20, 4)
        print('groundTs min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(groundTs), np.mean(groundTs), np.max(groundTs)))
        print('Percentage of negative gamma = ', 1e2*np.count_nonzero(gamma < 0)/n_tests)
        print('Percentage of negative groundTs = ', 1e2*np.count_nonzero(groundTs < 0)/n_tests)
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(15, 5)

    print('Squarings gain min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(squarings_gain),
                                                                np.mean(squarings_gain),
                                                                np.max(squarings_gain)))
    print('Percentage of negative squarings gain = ', 1e2*np.count_nonzero(squarings_gain < 0)/n_tests)
    print("#############################################################\n")

    fig.suptitle('{} Size {}, Number of test: {}'.format(title, N, n_tests))
    fig.canvas.set_window_title('{}-{}'.format(title, N))

    n, bins, patches = ax1.hist(gamma, 30, facecolor='g', alpha=0.75)
    ax1.set_xlabel('Gamma')
    ax1.set_ylabel('Frequency')
    ax1.title.set_text('Relative gain in L1 norm')
    ax1.grid(True)

    n, bins, patches = ax2.hist(squarings_gain, 20, facecolor='g', alpha=0.75)
    ax2.set_xlabel('Squarings gain')
    ax2.set_ylabel('Frequency')
    ax2.title.set_text('Number of squarings gained')
    ax2.grid(True)

    if (groundTs is not None):
        n, bins, patches = ax3.hist(groundTs, 30, facecolor='g', alpha=0.75)
        ax3.set_xlabel('Ground Truth')
        ax3.set_ylabel('Frequency')
        ax3.title.set_text('Difference from original norm')
        ax3.grid(True)

        n, bins, patches = ax4.hist(squarings_gain_GT, 30, facecolor='g', alpha=0.75)
        ax4.set_xlabel('Squarings gain')
        ax4.set_ylabel('Frequency')
        ax4.title.set_text('Number of squarings gained')
        ax4.grid(True)

    # plt.show()
    fig.savefig('scripts/plots/' + fig.canvas.get_window_title(), bbox_inches='tight')


whatToDo = True
forReal = True

# Section with tests on generated matrices
if whatToDo:
    np.random.seed(19680801)  # Fixing random state for reproducibility

    if forReal:
        MATRIX_SIZES = [4, 8, 16, 32, 64, 128, 256, 512, 1024]  # matrix size
        p = multiprocessing.Pool(psutil.cpu_count(logical=False))
        p.map(partial(run_random_tests), reversed(MATRIX_SIZES))
    else:
        run_random_tests(4)

# Section with tests on Toolbox matrices
else:
    matrices = importMatrices('scripts/testMatrices')
    analyseTheseMatrices(matrices)

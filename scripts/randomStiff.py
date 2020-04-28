import numpy as np
from numpy.linalg import norm
from testingStuff import testStuff
from usefulStuff import test_matrix, test_matrix_GT, maxnorm  # , generateStiffMatrix
from manageTestMatrices import james2014Generator, importMatrices
import multiprocessing
from functools import partial
from balanceMethods import new_balance, combinedOldNew
import psutil  # Used to get number of cores, physical instead of logical
from plotStuff import printData
from datetime import datetime

startTime = datetime.now()
np.set_printoptions(precision=2, linewidth=250, suppress=True)
testStuff()


def run_random_tests(matrix_size, n_tests=None):
    N = matrix_size
    # Do less test the bigger the matrix is
    if n_tests is None:
        n_tests = 102400 // N
    gamma = np.empty((n_tests, 2))
    groundTs = np.empty((n_tests, 2))
    squarings_gain = np.empty((n_tests, 2))
    squarings_gain_GT = np.empty((n_tests, 2))

    for k in range(0, n_tests):
        if(k % 500 == 0 and k != 0):
            print('Test {} for dimension {}'.format(k, N))

        while True:
            gt, A = james2014Generator(N)  # generateStiffMatrix(N)

            # Loop until to avoid creation of matrices that does not need balancing
            if norm(A, 1) > maxnorm:
                break

        # Be aware of what method is called here
        gamma[k][0], squarings_gain[k][0], groundTs[k][0], squarings_gain_GT[k][0] = test_matrix_GT(A, gt, new_balance)
        gamma[k][1], squarings_gain[k][1], groundTs[k][1], squarings_gain_GT[k][1] = test_matrix_GT(A, gt, combinedOldNew)

    # Saving results of multi-hour simulation
    np.savetxt('scripts/resultsData/gamma{}.tsv'.format(N), gamma, delimiter='\t')
    np.savetxt('scripts/resultsData/squarings_gain{}.tsv'.format(N), squarings_gain, delimiter='\t')
    np.savetxt('scripts/resultsData/groundTs{}.tsv'.format(N), groundTs, delimiter='\t')
    np.savetxt('scripts/resultsData/squarings_gain_GT{}.tsv'.format(N), squarings_gain_GT, delimiter='\t')

    # printData("NB", gamma[:, 0], squarings_gain[:, 0], groundTs[:, 0], squarings_gain_GT[:, 0], N)
    # printData("NBComb", gamma[:, 1], squarings_gain[:, 1], groundTs[:, 1], squarings_gain_GT[:, 1], N)


def analyseTheseMatrices(matrices):
    howMany = len(matrices)
    gamma = np.empty(howMany)
    squarings_gain = np.empty(howMany)

    for k in range(0, howMany):
        gamma[k], squarings_gain[k] = test_matrix(matrices[k])

    printData("Toolbox Matrices", gamma, squarings_gain, None, None, matrices[0].shape[0])

# ---------------------------------------------------------------------------------------------------------


whatToDo = True
forReal = True  # Around 3 hours and an half

# Section with tests on generated matrices
if whatToDo:
    np.random.seed(19680801)  # Fixing random state for reproducibility

    if forReal:
        MATRIX_SIZES = [4, 8, 16, 32, 64, 128, 256, 512, 1024]  # matrix size
        p = multiprocessing.Pool(psutil.cpu_count(logical=False))
        p.map(partial(run_random_tests), reversed(MATRIX_SIZES))
    else:
        run_random_tests(4, 1000)

# Section with tests on Toolbox matrices
else:
    matrices = importMatrices('scripts/testMatrices')
    analyseTheseMatrices(matrices)


print(datetime.now() - startTime)

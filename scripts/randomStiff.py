import numpy as np
from numpy.linalg import norm
from testingStuff import testStuff
from usefulStuff import test_matrix, test_matrix_GT, maxnorm  # , generateStiffMatrix
from manageTestMatrices import james2014Generator, importMatrices
import multiprocessing
from functools import partial
from balanceMethods import new_balance, combinedOldNew
import psutil  # Used to get number of cores, physical instead of logical
# from plotStuff import printData
from datetime import datetime
import os
import glob

startTime = datetime.now()
np.set_printoptions(precision=2, linewidth=250, suppress=True)
testStuff()


def run_random_tests(data, n_tests=None):
    print(data)
    N = data[0]
    np.random.seed(data[1])
    ID = os.getpid()
    print("Start " + str(ID))
    # Do less test the bigger the matrix is
    if n_tests is None:
        n_tests = (102400 // N) // psutil.cpu_count(logical=False)
    gamma = np.empty((n_tests, 2))
    groundTs = np.empty((n_tests, 2))
    squarings_gain = np.empty((n_tests, 2), dtype=np.int32)
    squarings_GT_gain = np.empty((n_tests, 2), dtype=np.int32)
    iterations = np.empty((n_tests, 2), dtype=np.int32)

    for k in range(0, n_tests):
        # if(k % (10240 // N) == 0 and k != 0):
        #     print('Test {} for dimension {}'.format(k, N))

        while True:
            gt, A = james2014Generator(N)  # generateStiffMatrix(N)

            # Loop until to avoid creation of matrices that does not need balancing
            if norm(A, 1) > maxnorm:
                break

        # Be aware of what method is called here
        # Upacking this to avoid a couple of long lines
        unpack = test_matrix_GT(A, gt, new_balance)

        gamma[k][0] = unpack[0]
        squarings_gain[k][0] = unpack[1]
        groundTs[k][0] = unpack[2]
        squarings_GT_gain[k][0] = unpack[3]
        iterations[k][0] = unpack[4]

        unpack = test_matrix_GT(A, gt, combinedOldNew)

        gamma[k][1] = unpack[0]
        squarings_gain[k][1] = unpack[1]
        groundTs[k][1] = unpack[2]
        squarings_GT_gain[k][1] = unpack[3]
        iterations[k][1] = unpack[4]

    # Saving results of multi-hour simulation
    np.savetxt('scripts/tempResults/gamma{}.tsv'.format(ID), gamma, delimiter='\t')
    np.savetxt('scripts/tempResults/squarings_gain{}.tsv'.format(ID), squarings_gain, delimiter='\t', fmt='%d')
    np.savetxt('scripts/tempResults/groundTs{}.tsv'.format(ID), groundTs, delimiter='\t')
    np.savetxt('scripts/tempResults/squarings_GT_gain{}.tsv'.format(ID), squarings_GT_gain, delimiter='\t', fmt='%d')
    np.savetxt('scripts/tempResults/iterations{}.tsv'.format(ID), iterations, delimiter='\t', fmt='%d')

    print("Stop " + str(ID))

    # printData("NB", gamma[:, 0], squarings_gain[:, 0], groundTs[:, 0], squarings_GT_gain[:, 0], N)
    # printData("NBComb", gamma[:, 1], squarings_gain[:, 1], groundTs[:, 1], squarings_GT_gain[:, 1], N)


def analyseTheseMatrices(matrices):
    howMany = len(matrices)
    gamma = np.empty(howMany)
    squarings_gain = np.empty(howMany)

    for k in range(0, howMany):
        gamma[k], squarings_gain[k] = test_matrix(matrices[k])

    # printData("Toolbox Matrices", gamma, squarings_gain, None, None, matrices[0].shape[0])


def emptyDirectory(d):
    filesToRemove = [os.path.join(d, f) for f in os.listdir(d)]
    for f in filesToRemove:
        os.remove(f)

# ---------------------------------------------------------------------------------------------------------


whatToDo = True

# Section with tests on generated matrices
if whatToDo:
    matrixSizes = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
    # matrixSizes = [4]
    resultsNames = [['gamma', np.float64, '%.18e'],
                    ['squarings_gain', np.int32, '%d'],
                    ['groundTs', np.float64, '%.18e'],
                    ['squarings_GT_gain', np.int32, '%d'],
                    ['iterations', np.int32, '%d']]
    nRealCores = psutil.cpu_count(logical=False)
    emptyDirectory(os.getcwd() + '/scripts/tempResults/')  # Start clean

    # Launch all cores on a given size
    for s in reversed(matrixSizes):
        arrS = [[s, i] for i in range(nRealCores)]  # Need to provide different seed to each process
        p = multiprocessing.Pool(nRealCores)
        p.map(partial(run_random_tests), arrS)
        p.close()
        p.join()

        # Read temp files and gather results
        for r in resultsNames:
            fileList = glob.glob('scripts/tempResults/{}*'.format(r[0]))
            merger = None
            for f in fileList:
                if merger is None:
                    merger = np.loadtxt(f, dtype=r[1])
                else:
                    merger = np.concatenate((merger, np.loadtxt(f, dtype=r[1])))

            np.savetxt('scripts/tempResultsDone/{}{}.tsv'.format(r[0], s), merger, delimiter='\t', fmt=r[2])

        emptyDirectory(os.getcwd() + '/scripts/tempResults/')


# Section with tests on Toolbox matrices
else:
    matrices = importMatrices('scripts/testMatrices')
    analyseTheseMatrices(matrices)


print(datetime.now() - startTime)

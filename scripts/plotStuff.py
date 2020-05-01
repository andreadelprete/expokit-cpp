import numpy as np
import matplotlib.pyplot as plt


def printData(title, gamma, squarings_gain, groundTs, squarings_GT_gain, N):
    n_tests = gamma.shape[0]
    # Console output section
    print('Matrix size = ', N)
    print('Gamma min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(gamma), np.mean(gamma), np.max(gamma)))
    if groundTs is not None and squarings_GT_gain is not None:
        fig, (ax1, ax3) = plt.subplots(1, 2)
        fig.set_size_inches(15, 5)
        print('groundTs min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(groundTs), np.mean(groundTs), np.max(groundTs)))
        print('Percentage of negative gamma = ', 1e2*np.count_nonzero(gamma < 0)/n_tests)
        print('Percentage of negative groundTs = ', 1e2*np.count_nonzero(groundTs < 0)/n_tests)
    else:
        fig, ax1 = plt.subplots(1, 1)
        fig.set_size_inches(15, 5)
    print('Squarings gain min-avg-max = %7.3f  %7.3f  %7.3f' % (np.min(squarings_gain),
                                                                np.mean(squarings_gain),
                                                                np.max(squarings_gain)))
    print('Percentage of negative squarings gain = ', 1e2*np.count_nonzero(squarings_gain < 0)/n_tests)
    print("#############################################################\n")

    # Column bar section
    fig.suptitle('{} Size {}, Number of test: {}'.format(title, N, n_tests))
    fig.canvas.set_window_title('{}-{}'.format(title, N))

    n, bins, patches = ax1.hist(gamma, 30, facecolor='g', alpha=0.75)
    ax1.set_xlabel('Gamma')
    ax1.set_ylabel('Frequency')
    ax1.title.set_text('Relative gain in L1 norm')
    ax1.grid(True)

    if (groundTs is not None):
        n, bins, patches = ax3.hist(groundTs, 30, facecolor='g', alpha=0.75)
        ax3.set_xlabel('Ground Truth')
        ax3.set_ylabel('Frequency')
        ax3.title.set_text('Difference from original norm')
        ax3.grid(True)

    fig.savefig('scripts/plots/' + fig.canvas.get_window_title(), bbox_inches='tight')


# Tables section
def createSquaringsTableColumn(squaringsGain):
    n_tests = squaringsGain.shape[0]
    result = np.empty(4, dtype='object')

    unique, counts = np.unique(squaringsGain, return_counts=True)
    d = dict(zip(unique.tolist(), counts.tolist()))

    for i in range(-1, 3):  # Make sure all options are considered
        if i not in d:
            result[i + 1] = "0.0%"
        else:
            result[i + 1] = "{:.1f}%".format(((d[i]/n_tests) * 100))

    return result


def createTable(matSizes):
    yLabels = np.array(['iter', '-1', '0', '+1', '+2'], dtype='object')
    tableShape = (5, len(matSizes))

    NBvsRodney = np.empty(tableShape, dtype='object')
    NBCombvsRodney = np.empty(tableShape, dtype='object')
    NBCombvsGT = np.empty(tableShape, dtype='object')

    for i in range(len(matSizes)):
        N = matSizes[i]
        squarings_gain = np.loadtxt('scripts/tempResultsDone/squarings_gain{}.tsv'.format(N), dtype=np.int32)
        squarings_GT_gain = np.loadtxt('scripts/tempResultsDone/squarings_GT_gain{}.tsv'.format(N), dtype=np.int32)
        iterations = np.loadtxt('scripts/tempResultsDone/iterations{}.tsv'.format(N), dtype=np.int32)

        NBvsRodney[0, i] = '{:.0f}'.format(np.mean(iterations[:, 0]).round())
        NBvsRodney[1:, i] = createSquaringsTableColumn(squarings_gain[:, 0])

        NBCombvsRodney[0, i] = '{:.0f}'.format(np.mean(iterations[:, 1]).round())
        NBCombvsRodney[1:, i] = createSquaringsTableColumn(squarings_gain[:, 1])

        NBCombvsGT[0, i] = NBCombvsRodney[0, i]
        NBCombvsGT[1:, i] = createSquaringsTableColumn(squarings_GT_gain[:, 1])

    listTable = ((NBvsRodney, "NBvsRodney"), (NBCombvsRodney, "NBCombvsRodney"), (NBCombvsGT, "NBCombvsGT"))

    for tableData in listTable:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_title(tableData[1])
        ax1.axis('off')
        ax1.table(tableData[0],
                  loc='center',
                  rowLabels=yLabels,
                  colLabels=matSizes)

        fig.savefig('scripts/tables/{}'.format(tableData[1]), bbox_inches='tight', dpi=600)


# ------------------------------------------------------------------------------
def plotEverything():
    matrixSizes = [4, 8, 16, 32, 64, 128, 256, 512, 1024]  # matrix size
    createTable(matrixSizes)

    for ms in matrixSizes:
        gamma = np.loadtxt('scripts/tempResultsDone/gamma{}.tsv'.format(ms))
        groundTs = np.loadtxt('scripts/tempResultsDone/groundTs{}.tsv'.format(ms))
        squarings_gain = np.loadtxt('scripts/tempResultsDone/squarings_gain{}.tsv'.format(ms), dtype=np.int32)
        squarings_GT_gain = np.loadtxt('scripts/tempResultsDone/squarings_GT_gain{}.tsv'.format(ms), dtype=np.int32)

        printData("NB", gamma[:, 0], squarings_gain[:, 0], groundTs[:, 0], squarings_GT_gain[:, 0], ms)
        printData("NBComb", gamma[:, 1], squarings_gain[:, 1], groundTs[:, 1], squarings_GT_gain[:, 1], ms)


plotEverything()

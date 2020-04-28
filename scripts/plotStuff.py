import numpy as np
import matplotlib.pyplot as plt


def printData(title, gamma, squarings_gain, groundTs, squarings_gain_GT, N):
    n_tests = gamma.shape[0]
    # Console output section
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

    # Column bar section
    fig.suptitle('{} Size {}, Number of test: {}'.format(title, N, n_tests))
    fig.canvas.set_window_title('{}-{}'.format(title, N))

    n, bins, patches = ax1.hist(gamma, 30, facecolor='g', alpha=0.75)
    ax1.set_xlabel('Gamma')
    ax1.set_ylabel('Frequency')
    ax1.title.set_text('Relative gain in L1 norm')
    ax1.grid(True)

    n, bins, patches = ax2.hist(squarings_gain, 20, facecolor='g', alpha=0.75)
    ax2.set_xlabel('Squarings gained')
    ax2.set_ylabel('Frequency')
    ax2.title.set_text('Squarings gain wrt Rodney balance')
    ax2.grid(True)

    if (groundTs is not None):
        n, bins, patches = ax3.hist(groundTs, 30, facecolor='g', alpha=0.75)
        ax3.set_xlabel('Ground Truth')
        ax3.set_ylabel('Frequency')
        ax3.title.set_text('Difference from original norm')
        ax3.grid(True)

        n, bins, patches = ax4.hist(squarings_gain_GT, 30, facecolor='g', alpha=0.75)
        ax4.set_xlabel('Squarings gained')
        ax4.set_ylabel('Frequency')
        ax4.title.set_text('Squarings gain wrt Ground Truth')
        ax4.grid(True)

    fig.savefig('scripts/plots/' + fig.canvas.get_window_title(), bbox_inches='tight')


# Tables section
def createTableColumn(squaringsGain, N):
    n_tests = squaringsGain.shape[0]
    result = np.empty(4, dtype='object')

    unique, counts = np.unique(squaringsGain, return_counts=True)
    d = dict(zip(unique.astype('int').tolist(), counts.tolist()))

    for i in range(-1, 3):  # Make sure all options are considered
        if i not in d:
            result[i + 1] = "0.0%"
        else:
            result[i + 1] = "{:.1f}%".format(((d[i]/n_tests) * 100))

    return result


def createTable(matSizes):
    yLabels = np.array(['-1', '0', '+1', '+2'], dtype='object')
    tableShape = (4, len(matSizes))

    NBvsRodney = np.empty(tableShape, dtype='object')
    NBCombvsRodney = np.empty(tableShape, dtype='object')
    NBCombvsGT = np.empty(tableShape, dtype='object')

    for i in range(len(matSizes)):
        N = matSizes[i]
        squarings_gain = np.loadtxt('scripts/resultsData/squarings_gain{}.tsv'.format(N))
        squarings_gain_GT = np.loadtxt('scripts/resultsData/squarings_gain_GT{}.tsv'.format(N))

        NBvsRodney[:, i] = createTableColumn(squarings_gain[:, 0], N)
        NBCombvsRodney[:, i] = createTableColumn(squarings_gain[:, 1], N)
        NBCombvsGT[:, i] = createTableColumn(squarings_gain_GT[:, 1], N)

    listTable = ((NBvsRodney, "NBvsRodney"), (NBCombvsRodney, "NBCombvsRodney"), (NBCombvsGT, "NBCombvsGT"))

    for tableData in listTable:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(tableData[1])
        ax.axis('off')
        ax.table(tableData[0],
                 loc='center',
                 rowLabels=yLabels,
                 colLabels=matSizes)
        # table.set_fontsize(10)
        # table.scale(1, 4)
        # fig.set_size_inches(35, 15)
        # fig.suptitle(tableData[1])
        fig.savefig('scripts/tables/{}'.format(tableData[1]), bbox_inches='tight', dpi=600)


# ------------------------------------------------------------------------------
'''
gamma = np.loadtxt('scripts/resultsData/gamma4.tsv')
groundTs = np.loadtxt('scripts/resultsData/groundTs4.tsv')
squarings_gain = np.loadtxt('scripts/resultsData/squarings_gain4.tsv')
squarings_gain_GT = np.loadtxt('scripts/resultsData/squarings_gain_GT4.tsv')

N = 4
# printData("NB", gamma[:, 0], squarings_gain[:, 0], groundTs[:, 0], squarings_gain_GT[:, 0], N)
# printData("NBComb", gamma[:, 1], squarings_gain[:, 1], groundTs[:, 1], squarings_gain_GT[:, 1], N)
'''

MATRIX_SIZES = [4, 8, 16, 32, 64, 128, 256, 512, 1024]  # matrix size
createTable(MATRIX_SIZES)

#ifndef BALANCING_METHODS_H
#define BALANCING_METHODS_H

#include <Eigen/Core>
#include "utils/stop-watch.h"

using namespace Eigen;


namespace expokit {

template <typename T, int N>
class BalancingMethods {
    int size;

    // Typedefs to make code more readable
    typedef Matrix<T, N, 1> VectorType;
    typedef Matrix<T, N, N> MatrixType;

    typedef const Ref<const VectorType> RefVector;
    typedef const Ref<const MatrixType> RefMatrix;

    typedef Ref<MatrixType> RefOutMatrix;
    typedef Ref<VectorType> RefOutVector;

    Matrix<T, 1, N> columnNorms; // Used in balancing
    MatrixType allColumnNorms; // Used in balancing
    MatrixType temp; // Used in combined balancing

    void init(int n);

public:
    BalancingMethods();
    explicit BalancingMethods(int n);

    int balanceNew(RefMatrix& A, RefOutMatrix B, RefOutMatrix D, RefOutMatrix Dinv, int maxIter = 0);
    int balanceRodney(RefMatrix& A, RefOutMatrix B, RefOutMatrix D, RefOutMatrix Dinv, int maxIter = 0);
    int balanceCombined(RefMatrix& A, RefOutMatrix B, RefOutMatrix D, RefOutMatrix Dinv);

};

template <typename T, int N>
BalancingMethods<T, N>::BalancingMethods()
{
    if (N == Dynamic) {
        init(2); // Init to dym 2 for no particular reason
    } else {
        init(N);
    }
}

template <typename T, int N>
BalancingMethods<T, N>::BalancingMethods(int n)
{
    init(n);
}

template <typename T, int N>
void BalancingMethods<T, N>::init(int n)
{
    size = n;
    columnNorms.resize(1, size);
    allColumnNorms = MatrixType::Zero(size, size);
    temp.resize(n, n);
}

template <typename T, int N>
int BalancingMethods<T, N>::balanceNew(RefMatrix& A, RefOutMatrix B, RefOutMatrix D, RefOutMatrix Dinv, int maxIter)
{
    // START_PROFILER("preamble");
    B = A; // Copy
    D = MatrixType::Identity(size, size);
    Dinv = MatrixType::Identity(size, size);

    columnNorms = B.cwiseAbs().colwise().sum();
    allColumnNorms = MatrixType::Zero(size, size);
    int jMax = 0;
    columnNorms.maxCoeff(&jMax);

    T vMin = columnNorms(jMax); // Do not do worse than this
    int iMin = 0;

    if (maxIter == 0)
        maxIter = 15 * size; // Simplyfied - to be checked

    // STOP_PROFILER("preamble");

    int it = 0;
    for (; it < maxIter; ++it) {
        // START_PROFILER("innerfor");
        for (int k = 0; k < size; ++k) {
            T possibleNorm;
            if (k == jMax) {
                allColumnNorms.row(k) = columnNorms + B.row(k).cwiseAbs();
                allColumnNorms(k, k) = (columnNorms(k) + abs(B(k, k))) / 2;
                possibleNorm = allColumnNorms.row(k).maxCoeff();
            } else {
                allColumnNorms.row(k) = columnNorms - (B.row(k).cwiseAbs() / 2);
                allColumnNorms(k, k) = (columnNorms(k) * 2) - abs(B(k, k));
                possibleNorm = allColumnNorms.row(k).maxCoeff();
            }

            if (possibleNorm < vMin) {
                vMin = possibleNorm;
                iMin = k;
            }
        }
        // STOP_PROFILER("innerfor");

        // START_PROFILER("confirm");
        // If last loop failed to find better solution
        if (vMin >= columnNorms(jMax))
            break; // Do not continue looking

        if (iMin == jMax) {
            D(iMin, iMin) /= 2;
            Dinv(iMin, iMin) *= 2;
            B.col(iMin) /= 2;
            B.row(iMin) *= 2;
        } else {
            D(iMin, iMin) *= 2;
            Dinv(iMin, iMin) /= 2;
            B.col(iMin) *= 2;
            B.row(iMin) /= 2;
        }

        // Save selected norms for new loop
        columnNorms = allColumnNorms.row(iMin);
        columnNorms.maxCoeff(&jMax);
        // STOP_PROFILER("confirm");
    }

    return it;
}

template <typename T, int N>
int BalancingMethods<T, N>::balanceRodney(RefMatrix& A, RefOutMatrix B, RefOutMatrix D, RefOutMatrix Dinv, int maxIter)
{
    D = MatrixType::Identity(size, size);
    Dinv = MatrixType::Identity(size, size);

    B = A;

    T c, r, s, f;

    bool converged = false;
    int it = 0;
    while (!converged) {
        converged = true;

        for (int i = 0; i < size; ++i) {
            c = B.col(i).template lpNorm<1>();
            r = B.row(i).template lpNorm<1>();
            s = c * c + r * r;
            f = 1;

            while (c < r / 2 && c != 0) {
                c = c * 2;
                r = r / 2;
                f = f * 2;
            }

            while (c >= r * 2 && r != 0) {
                c = c / 2;
                r = r * 2;
                f = f / 2;
            }

            if (c * c + r * r < 0.95 * s) {
                converged = false;
                D(i, i) = f * D(i, i);
                Dinv(i, i) = Dinv(i, i) / f;
                B.col(i) = f * B.col(i);
                B.row(i) = B.row(i) / f;
            }
        }
        it++;
        if (maxIter != 0 && it > maxIter)
            break;
    }

    return it;
}



template <typename T, int N>
int BalancingMethods<T, N>::balanceCombined(RefMatrix& A, RefOutMatrix B, RefOutMatrix D, RefOutMatrix Dinv) {
    T oNorm = A.cwiseAbs().colwise().sum().maxCoeff();
    int it = balanceRodney(A, temp, D, Dinv);
    T normRodney = temp.cwiseAbs().colwise().sum().maxCoeff();

    // Compute again from scratch
    if(oNorm < normRodney) {
        it = balanceNew(A, B, D, Dinv);
    } else { // chain results
        it += balanceNew(temp, B, D, Dinv);
    }

    return it;
}


}

#endif
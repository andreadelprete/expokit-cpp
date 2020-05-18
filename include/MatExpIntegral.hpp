#ifndef MATEXPINTEGRAL_H
#define MATEXPINTEGRAL_H

#include "MatrixExponential.hpp"
#include <Eigen/Core>

using namespace Eigen;

namespace expokit {

template <typename T>
class MatExpIntegral {

private:
    int size;
    int size2;

    // Typedefs to make code more readable
    typedef Matrix<T, Dynamic, 1> VectorType;
    typedef Matrix<T, Dynamic, Dynamic> MatrixType;

    typedef const Ref<const VectorType> RefVector;
    typedef const Ref<const MatrixType> RefMatrix;

    typedef Ref<MatrixType> RefOutMatrix;
    typedef Ref<VectorType> RefOutVector;

    MatrixExponential<T, Dynamic> expUtil1;
    MatrixExponential<T, Dynamic> expUtil2;

    MatrixType C, zInt, temp1, temp2, temp3, temp4, ex;

    void init(int n);

public:
    void computeExpIntegral(RefMatrix& A, RefOutMatrix out, T t, T dt = 0);

    MatExpIntegral(int n);
    MatExpIntegral();
};

template <typename T>
MatExpIntegral<T>::MatExpIntegral()
{
    // if (N == Dynamic) {
        init(2); // Init to dym 2 for no particular reason
    // } else {
    //     init(N);
    // }
}

template <typename T>
MatExpIntegral<T>::MatExpIntegral(int n)
{
    init(n);
}

template <typename T>
void MatExpIntegral<T>::init(int n)
{
    size = n;
    size2 = n * 2;
    expUtil1.resize(size);
    expUtil2.resize(size2);

    C = MatrixType::Zero(size2, size2);
    C.block(0, size, size, size) = MatrixType::Identity(size, size);

    zInt = MatrixType::Zero(size2, size);
    zInt.block(size, 0, size, size) = MatrixType::Identity(size, size);

    temp1.resize(size2, size2);

    temp2.resize(size2, size);

    temp3.resize(size, size);

    temp4.resize(size2, size2);

    ex.resize(size, size);
}

template <typename T>
void MatExpIntegral<T>::computeExpIntegral(RefMatrix& A, RefOutMatrix out, T t, T dt)
{
    if (dt != 0) {
        int n = int(t / dt);

        out = MatrixType::Zero(size, size);
        for (int i = 0; i < n; ++i) {
            temp3 = i * dt * A;
            expUtil1.compute(temp3, ex);
            out += ex;
        }

        return;
    }

    C.block(0, 0, size, size) = A;

    temp4 = t * C;
    expUtil2.compute(temp4, temp1);

    temp2.noalias() = temp1 * zInt;

    out = temp2.block(0, 0, size, size);
}

}
#endif
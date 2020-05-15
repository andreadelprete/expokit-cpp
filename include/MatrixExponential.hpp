/*
    Insert fancy disclaimer here
*/

#ifndef EIGEN_MATRIX_EXPONENTIAL_EXPOKIT
#define EIGEN_MATRIX_EXPONENTIAL_EXPOKIT

#ifdef EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
#define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
#define EIGEN_MALLOC_ALLOWED
#define EIGEN_MALLOC_NOT_ALLOWED
#endif

#include "unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"
#include "utils/stop-watch.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
// you need to define MAIN__ function because f2c complains if it doesn't find one
int MAIN__() { return 0; }
}
#endif

#define DEFAULT_MIN_SQUARINGS 10
#define METAPROD_SIZE 50

using namespace Eigen;

namespace expokit {

template <typename T, int N>
class MatrixExponential {

private:
    int size;

    // Typedefs to make code more readable
    typedef Matrix<T, N, 1> VectorType;
    typedef Matrix<T, N, N> MatrixType;

    typedef const Ref<const VectorType> RefVector;
    typedef const Ref<const MatrixType> RefMatrix;

    typedef Ref<MatrixType> RefOutMatrix;
    typedef Ref<VectorType> RefOutVector;

    MatrixType U, V, numer, denom, A_scaled, A2, A4, A6, A8, tmp, eye, tmp2, D, Dinv, Abal;
    VectorType v_tmp, vTmp1, vTmp2;
    PartialPivLU<MatrixType> ppLU;
    int squarings;

    MatrixType metaProds[METAPROD_SIZE]; // Should be enough to not check every time?
    MatrixType prevA, deltaA; // Used in delta update
    bool delta, deltaUsed;
    int minSquarings; // Default is 10

public:
    MatrixExponential();
    explicit MatrixExponential(int n);

    void resize(int n);
    void resetDelta();

    int getSquarings() const { return squarings; }

    void useDelta(bool delta) { this->delta = delta; }
    bool isDeltaUsed() { return delta; }
    bool wasDeltaUsed() { return deltaUsed; }

    int getMinSquarings() { return minSquarings; }
    void setMinSquarings(int minSquarings) { this->minSquarings = minSquarings; }

    /** 
     * Compute the exponential of the given matrix arg and writes it in result.
     */
    void compute(RefMatrix& A, RefOutMatrix out);
    // void balanceCompute(RefMatrix& A, RefOutMatrix out);

    /** Compute the product between the exponential of the given matrix arg and the given
     * vector v. The result is written it the output variable result.
     * The optional parameter vec_squarings specifies how many of the squaring operations
     * are performed through matrix-vector products. The remaining squaring operations are
     * then performed through matrix-matrix products, as in the classical scaling-and-squaring
     * algorithm. Therefore, specifying vec_squarings=0 this method behaves in the same way
     * as the compute method. Specifying vec_squarings<0, the number of vec_squarings is
     * automatically computed, but it may not be the best one, so if the user wants to achieve
     * maximum speed, she/he should test different values of vec_squarings.
     */
    void computeExpTimesVector(RefMatrix& A, RefVector& v, RefOutVector out, int vec_squarings = -1);
    // void balanceComputeExpTimesVector(RefMatrix& A, RefVector& v, RefOutVector out, int vec_squarings);
    //void computeExpTimesVector(RefVector v, , int vec_squarings = -1);


private:
    void init(int n);


    int determineSquarings(const double l1norm);

    void computeUV(RefMatrix& A);
    // Overloading that does not compute squarings inside
    void computeUV(RefMatrix& A, int squarings);

    /** \brief Compute the (3,3)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade3(RefMatrix& A);

    /** \brief Compute the (5,5)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade5(RefMatrix& A);

    /** \brief Compute the (7,7)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade7(RefMatrix& A);

    /** \brief Compute the (9,9)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade9(RefMatrix& A);

    /** \brief Compute the (13,13)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade13(RefMatrix& A);
};

template <typename T, int N>
MatrixExponential<T, N>::MatrixExponential()
{
    if (N == Dynamic) {
        init(2); // Init to dym 2 for no particular reason
    } else {
        init(N);
    }
}

template <typename T, int N>
MatrixExponential<T, N>::MatrixExponential(int n)
{
    init(n);
}

template <typename T, int N>
void MatrixExponential<T, N>::init(int n)
{
    size = n;
    U.resize(n, n);
    V.resize(n, n);
    numer.resize(n, n);
    denom.resize(n, n);
    A2.resize(n, n);
    A4.resize(n, n);
    A6.resize(n, n);
    A8.resize(n, n);
    tmp.resize(n, n);
    tmp2.resize(n, n);
    A_scaled.resize(n, n);
    D.resize(n, n);
    Dinv.resize(n, n);
    Abal.resize(n, n);
    eye = MatrixType::Identity(n, n);
    ppLU = PartialPivLU<MatrixType>(n);
    v_tmp.resize(n);
    vTmp1.resize(n);
    vTmp2.resize(n);
    squarings = 0;
    minSquarings = DEFAULT_MIN_SQUARINGS;
    delta = false;
    for (int i = 0; i < METAPROD_SIZE; i++) {
        metaProds[i].resize(n, n);
    }
    resetDelta();


}

template <typename T, int N>
void MatrixExponential<T, N>::resize(int n)
{
    init(n);
    resetDelta();
}

template <typename T, int N>
void MatrixExponential<T, N>::resetDelta()
{
    // Necessary to trigger full computation at the first call
    prevA = MatrixType::Zero(size, size);
}

template <typename T, int N>
void MatrixExponential<T, N>::compute(RefMatrix& A, RefOutMatrix out)
{
    deltaUsed = false;
    int deltaSquarings;
    if (delta) {
        deltaA = A - prevA;
        // TODO think about directly the method for the 1 norm
        const double l1normDeltaA = deltaA.cwiseAbs().colwise().sum().maxCoeff();
        // Speedup only if the difference in number of squaring
        // from the complete recomputation is greater than 2
        deltaSquarings = determineSquarings(l1normDeltaA);

        deltaUsed = (squarings >= deltaSquarings);
    }

    if (deltaUsed) {
        computeUV(deltaA, deltaSquarings);
        numer = U + V;
        denom = -U + V;
        ppLU.compute(denom);
        // There wuold be alias with squarings = 0, but it's not possible by logic
        metaProds[squarings] = metaProds[squarings] * ppLU.solve(numer);
    } else {
        computeUV(A);
        numer = U + V;
        denom = -U + V;
        ppLU.compute(denom);
        metaProds[squarings] = ppLU.solve(numer);
    }

    // undo scaling by repeated squaring
    for (int i = squarings; i > 0; --i) {
        metaProds[i - 1].noalias() = metaProds[i] * metaProds[i];
    }

    out = metaProds[0];
    prevA = A;
}

/*
template <typename T, int N>
void MatrixExponential<T, N>::balanceCompute(RefMatrix& A, RefOutMatrix out)
{
    balanceRodney(A);
    // std::cout << Abal << std::endl;
    computeUV(Abal);
    numer = U + V;
    denom = -U + V;
    ppLU.compute(denom);
    metaProds[squarings] = ppLU.solve(numer);

    for (int i = squarings; i > 0; --i) {
        metaProds[i - 1].noalias() = metaProds[i] * metaProds[i];
    }

    tmp.noalias() = D * metaProds[0];
    out.noalias() = tmp * Dinv;
}*/

template <typename T, int N>
void MatrixExponential<T, N>::computeExpTimesVector(RefMatrix& A, RefVector& v, RefOutVector out, int vec_squarings)
{
    computeUV(A); // Pade approximant is (U+V) / (-U+V)
    numer = U + V;
    denom = -U + V;
    ppLU.compute(denom);
    tmp = ppLU.solve(numer);

    // undo scaling by repeated squaring
    /*
      * If the matrix size is n and the number of squaring is s, applying
      * matrix-vector multiplications is only convenient if
      *             n > (2^s)/s
      *
      * This is (2^s)/s for s in range(1,10): array([2, 2, 2, 4, 6, 10, 18, 32, 56])
      */
    if (vec_squarings < 0) {
        int n = (int)v.size();
        const int b[] = { 2, 2, 2, 4, 6, 10, 18, 32, 56 };
        for (int i = 8; i >= 0; i--)
            if (n > b[i])
                vec_squarings = i + 1;
    }

    // number of squarings implemented via matrix-matrix multiplications
    int mat_squarings = squarings - vec_squarings;
    if (mat_squarings < 0) {
        mat_squarings = 0;
        vec_squarings = squarings;
    }

    for (int i = 0; i < mat_squarings; i++) {
        tmp2.noalias() = tmp * tmp;
        tmp = tmp2;
    }

    //int two_pow_s = (int)std::pow(2, vec_squarings);
    unsigned int two_pow_s = 1U << (unsigned int)vec_squarings;

    v_tmp = v;
    for (unsigned int i = 0; i < two_pow_s; i++) {
        out.noalias() = tmp * v_tmp;
        v_tmp = out;
    }
}

/*
template <typename T, int N>
void MatrixExponential<T, N>::balanceComputeExpTimesVector(RefMatrix& A, RefVector& v, RefOutVector out, int vec_squarings)
{
    balanceRodney(A);
    vTmp1.noalias() = Dinv * v;
    // std::cout << v << std::endl;
    computeExpTimesVector(Abal, vTmp1, vTmp2, vec_squarings);
    out.noalias() = D * vTmp2;
}*/

template <typename T, int N>
int MatrixExponential<T, N>::determineSquarings(const double l1norm)
{
    const double maxnorm = 5.371920351148152;
    int squars;
    std::frexp(l1norm / maxnorm, &squars);
    if (squars < 0)
        squars = 0;

    if (delta && squars < minSquarings)
        squars = minSquarings;

    return squars;
}

template <typename T, int N>
void MatrixExponential<T, N>::computeUV(RefMatrix& A)
{
    const double l1norm = A.cwiseAbs().colwise().sum().maxCoeff();
    squarings = 0;
    if (l1norm > 2.097847961257068e+000) {
        squarings = determineSquarings(l1norm);
        A_scaled = A.unaryExpr(Eigen::internal::MatrixExponentialScalingOp<double>(squarings));
        matrix_exp_pade13(A_scaled);
    } else if (l1norm < 1.495585217958292e-002) {
        matrix_exp_pade3(A);
    } else if (l1norm < 2.539398330063230e-001) {
        matrix_exp_pade5(A);
    } else if (l1norm < 9.504178996162932e-001) {
        matrix_exp_pade7(A);
    } else {
        matrix_exp_pade9(A);
    }
}

template <typename T, int N>
void MatrixExponential<T, N>::computeUV(RefMatrix& A, int squarings)
{
    // Pade approximant is (U+V) / (-U+V)
    this->squarings = squarings;
    A_scaled = A.unaryExpr(Eigen::internal::MatrixExponentialScalingOp<double>(squarings));
    const double l1norm = A_scaled.cwiseAbs().colwise().sum().maxCoeff();
    if (l1norm > 2.097847961257068e+000) {
        matrix_exp_pade13(A_scaled);
    } else if (l1norm < 1.495585217958292e-002) {
        matrix_exp_pade3(A_scaled);
    } else if (l1norm < 2.539398330063230e-001) {
        matrix_exp_pade5(A_scaled);
    } else if (l1norm < 9.504178996162932e-001) {
        matrix_exp_pade7(A_scaled);
    } else {
        matrix_exp_pade9(A_scaled);
    }
}

template <typename T, int N>
void MatrixExponential<T, N>::matrix_exp_pade3(RefMatrix& A)
{
    //    typedef typename Eigen::internal::NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 120., 60., 12., 1. };
    A2.noalias() = A * A;
    tmp.noalias() = b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V.noalias() = b[2] * A2 + b[0] * eye;
}

template <typename T, int N>
void MatrixExponential<T, N>::matrix_exp_pade5(RefMatrix& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 30240., 15120., 3360., 420., 30., 1. };
    A2.noalias() = A * A;
    A4.noalias() = A2 * A2;
    tmp = b[5] * A4 + b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V = b[4] * A4 + b[2] * A2 + b[0] * eye;
}

template <typename T, int N>
void MatrixExponential<T, N>::matrix_exp_pade7(RefMatrix& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1. };
    A2.noalias() = A * A;
    A4.noalias() = A2 * A2;
    A6.noalias() = A4 * A2;
    tmp = b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V = b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * eye;
}

template <typename T, int N>
void MatrixExponential<T, N>::matrix_exp_pade9(RefMatrix& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 17643225600., 8821612800., 2075673600., 302702400., 30270240.,
        2162160., 110880., 3960., 90., 1. };
    A2.noalias() = A * A;
    A4.noalias() = A2 * A2;
    A6.noalias() = A4 * A2;
    A8.noalias() = A6 * A2;
    tmp = b[9] * A8 + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V = b[8] * A8 + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * eye;
}

template <typename T, int N>
void MatrixExponential<T, N>::matrix_exp_pade13(RefMatrix& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 64764752532480000., 32382376266240000., 7771770303897600.,
        1187353796428800., 129060195264000., 10559470521600., 670442572800.,
        33522128640., 1323241920., 40840800., 960960., 16380., 182., 1. };

    A2.noalias() = A * A;
    A4.noalias() = A2 * A2;
    A6.noalias() = A4 * A2;
    V = b[13] * A6 + b[11] * A4 + b[9] * A2; // used for temporary storage
    tmp.noalias() = A6 * V;
    tmp += b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * eye; // this is a sum so no alias by default
    U.noalias() = A * tmp;
    tmp = b[12] * A6 + b[10] * A4 + b[8] * A2;
    V.noalias() = A6 * tmp;
    V += b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * eye;
}


} // end namespace Eigen

#endif // EIGEN_MATRIX_EXPONENTIAL

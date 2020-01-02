/*
    Insert fancy disclaimer here
*/

#ifndef EIGEN_MATRIX_EXPONENTIAL_EXPOKIT
#define EIGEN_MATRIX_EXPONENTIAL_EXPOKIT

#include "utils/stop-watch.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <stdio.h>
#include "unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"


using namespace Eigen;

namespace expokit {

template <typename type, int N>
class MatrixExponential {

private:
    // Typedefs to make code more readable
    typedef const Ref<const Matrix<type, N, 1> > RefVector;
    typedef const Ref<const Matrix<type, N, N> > RefMatrix;
    typedef Matrix<type, N, 1> VectorType;
    typedef Matrix<type, N, N> MatrixType;

    MatrixType U, V, numer, denom;
    MatrixType A_scaled, A2, A4, A6, A8, tmp, eye, tmp2;
    VectorType v_tmp;
    PartialPivLU<MatrixType> ppLU;
    int squarings;

public:
    MatrixExponential();

    int get_squarings() const
    {
        return squarings;
    }

    /** Compute the exponential of the given matrix arg and writes it in result.
     */
    Matrix<type, N, N> compute(const RefMatrix& arg);

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
    Matrix<type, N, 1> computeExpTimesVector(const RefMatrix& arg, const RefVector& v, int vec_squarings = -1);

private:
    void computeUV(const MatrixType& arg);

    /** \brief Compute the (3,3)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade3(const MatrixType& A);
    
    /** \brief Compute the (5,5)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade5(const MatrixType& A);

    /** \brief Compute the (7,7)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade7(const MatrixType& A);

    /** \brief Compute the (9,9)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade9(const MatrixType& A);

    /** \brief Compute the (13,13)-Pad&eacute; approximant to the exponential.
    *
    *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
    *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
    */
    void matrix_exp_pade13(const MatrixType& A);
};

template <typename type, int N>
MatrixExponential<type, N>::MatrixExponential()
{
    START_PROFILER("MatrixExponential::constructor");
    eye = MatrixType::Identity(N, N);
    ppLU = PartialPivLU<MatrixType>(N);
    squarings = 0;
    STOP_PROFILER("MatrixExponential::constructor");
}

template <typename type, int N>
Matrix<type, N, N> MatrixExponential<type, N>::compute(const RefMatrix& arg)
{
    START_PROFILER("MatrixExponential::compute");
    START_PROFILER("MatrixExponential::computeUV");
    computeUV(arg); // Pade approximant is (U+V) / (-U+V)
    STOP_PROFILER("MatrixExponential::computeUV");
    numer = U + V;
    denom = -U + V;
    START_PROFILER("MatrixExponential::computeLUdecomposition");
    ppLU.compute(denom);
    STOP_PROFILER("MatrixExponential::computeLUdecomposition");
    START_PROFILER("MatrixExponential::LUsolve");
    tmp = ppLU.solve(numer);
    STOP_PROFILER("MatrixExponential::LUsolve");
    //    std::cout<<"Squaring: "<<squarings<<std::endl;

    // undo scaling by repeated squaring
    START_PROFILER("MatrixExponential::squaring");
    for (int i = 0; i < squarings; i++) {
        tmp2.noalias() = tmp * tmp;
        tmp = tmp2;
    }
    STOP_PROFILER("MatrixExponential::squaring");
    return tmp;
    STOP_PROFILER("MatrixExponential::compute");
}

template <typename type, int N>
Matrix<type, N, 1> MatrixExponential<type, N>::computeExpTimesVector(const RefMatrix& arg, const RefVector& v, int vec_squarings)
{
    START_PROFILER("MatrixExponential::computeExpTimesVector");
    START_PROFILER("MatrixExponential:computeExpTimesVector:computeUV");
    computeUV(arg); // Pade approximant is (U+V) / (-U+V)
    STOP_PROFILER("MatrixExponential:computeExpTimesVector:computeUV");
    numer = U + V;
    denom = -U + V;
    START_PROFILER("MatrixExponential:computeExpTimesVector:computeLUdecomposition");
    ppLU.compute(denom);
    STOP_PROFILER("MatrixExponential:computeExpTimesVector:computeLUdecomposition");
    START_PROFILER("MatrixExponential:computeExpTimesVector:LUsolve");
    tmp = ppLU.solve(numer);
    STOP_PROFILER("MatrixExponential:computeExpTimesVector:LUsolve");
    //    std::cout<<"Squaring: "<<squarings<<std::endl;

    // undo scaling by repeated squaring
    START_PROFILER("MatrixExponential:computeExpTimesVector:squaringVector");
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
    int two_pow_s = (unsigned int) 1 << vec_squarings;

    v_tmp = v;
    Matrix<type, N, 1> result;
    for (int i = 0; i < two_pow_s; i++) {
        result.noalias() = tmp * v_tmp;
        v_tmp = result;
    }
    
    return result; // Relying on NRVO
    STOP_PROFILER("MatrixExponential:computeExpTimesVector:squaringVector");
    STOP_PROFILER("MatrixExponential::computeExpTimesVector");
}

template <typename type, int N>
void MatrixExponential<type, N>::computeUV(const MatrixType& arg)
{
    const double l1norm = arg.cwiseAbs().colwise().sum().maxCoeff();
    squarings = 0;
    if (l1norm < 1.495585217958292e-002) {
        matrix_exp_pade3(arg);
    } else if (l1norm < 2.539398330063230e-001) {
        matrix_exp_pade5(arg);
    } else if (l1norm < 9.504178996162932e-001) {
        matrix_exp_pade7(arg);
    } else if (l1norm < 2.097847961257068e+000) {
        matrix_exp_pade9(arg);
    } else {
        const double maxnorm = 5.371920351148152;
        std::frexp(l1norm / maxnorm, &squarings);
        if (squarings < 0)
            squarings = 0;
        A_scaled = arg.unaryExpr(Eigen::internal::MatrixExponentialScalingOp<double>(squarings));
        START_PROFILER("MatrixExponential::matrix_exp_pade13");
        matrix_exp_pade13(A_scaled);
        STOP_PROFILER("MatrixExponential::matrix_exp_pade13");
    }
}

template <typename type, int N>
void MatrixExponential<type, N>::matrix_exp_pade3(const MatrixType& A)
{
    //    typedef typename Eigen::internal::NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 120., 60., 12., 1. };
    A2.noalias() = A * A;
    tmp.noalias() = b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V.noalias() = b[2] * A2 + b[0] * eye;
}

template <typename type, int N>
void MatrixExponential<type, N>::matrix_exp_pade5(const MatrixType& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 30240., 15120., 3360., 420., 30., 1. };
    A2 = A * A;
    A4 = A2 * A2;
    tmp = b[5] * A4 + b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V = b[4] * A4 + b[2] * A2 + b[0] * eye;
}

template <typename type, int N>
void MatrixExponential<type, N>::matrix_exp_pade7(const MatrixType& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1. };
    A2 = A * A;
    A4 = A2 * A2;
    A6 = A4 * A2;
    tmp = b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V = b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * eye;
}

template <typename type, int N>
void MatrixExponential<type, N>:: matrix_exp_pade9(const MatrixType& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 17643225600., 8821612800., 2075673600., 302702400., 30270240.,
        2162160., 110880., 3960., 90., 1. };
    A2 = A * A;
    A4 = A2 * A2;
    A6 = A4 * A2;
    A8 = A6 * A2;
    tmp = b[9] * A8 + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * eye;
    U.noalias() = A * tmp;
    V = b[8] * A8 + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * eye;
}

template <typename type, int N>
void MatrixExponential<type, N>:: matrix_exp_pade13(const MatrixType& A)
{
    //    typedef typename NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
    typedef double RealScalar;
    const RealScalar b[] = { 64764752532480000., 32382376266240000., 7771770303897600.,
        1187353796428800., 129060195264000., 10559470521600., 670442572800.,
        33522128640., 1323241920., 40840800., 960960., 16380., 182., 1. };
    A2.noalias() = A * A;
    A4.noalias() = A2 * A2;
    A6.noalias() = A4 * A2;
    V.noalias() = b[13] * A6;
    V.noalias() += b[11] * A4;
    V.noalias() += b[9] * A2; // used for temporary storage
    tmp.noalias() = A6 * V;
    tmp.noalias() += b[7] * A6;
    tmp.noalias() += b[5] * A4;
    tmp.noalias() += b[3] * A2;
    tmp.noalias() += b[1] * eye;
    U.noalias() = A * tmp;
    tmp.noalias() = b[12] * A6;
    tmp.noalias() += b[10] * A4;
    tmp.noalias() += b[8] * A2;
    V.noalias() = A6 * tmp;
    V.noalias() += b[6] * A6;
    V.noalias() += b[4] * A4;
    V.noalias() += b[2] * A2;
    V.noalias() += b[0] * eye;
}

} // end namespace Eigen

#endif // EIGEN_MATRIX_EXPONENTIAL

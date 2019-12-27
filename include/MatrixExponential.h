// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009, 2010, 2013 Jitse Niesen <jitse@maths.leeds.ac.uk>
// Copyright (C) 2011, 2013 Chen-Pang He <jdh8@ms63.hinet.net>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_MATRIX_EXPONENTIAL_EXPOKIT
#define EIGEN_MATRIX_EXPONENTIAL_EXPOKIT

//#include "StemFunction.h"
#include "unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"
#include "unsupported/Eigen/src/MatrixFunctions/StemFunction.h"
#include <iostream>
#include <stdio.h>

#include "utils/stop-watch.h"

namespace Eigen {

template <typename MatrixType, typename VectorType>
class MatrixExponential {
    MatrixType U, V, numer, denom;
    MatrixType A_scaled, A2, A4, A6, A8, tmp, eye, tmp2;
    VectorType v_tmp;
    PartialPivLU<MatrixType> ppLU;
    int squarings;

public:
    MatrixExponential(int n)
    {
        START_PROFILER("MatrixExponential::resize");
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
        eye = MatrixType::Identity(n, n);
        ppLU = PartialPivLU<MatrixType>(n);
        v_tmp.resize(n);
        STOP_PROFILER("MatrixExponential::resize");
    }

    int get_squarings() const
    {
        return squarings;
    }

    void compute(const MatrixType& arg, MatrixType& result)
    {
        START_PROFILER("MatrixExponential::compute");
        START_PROFILER("MatrixExponential::computeUV");
        computeUV(arg, U, V, squarings); // Pade approximant is (U+V) / (-U+V)
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
        result = tmp;
        STOP_PROFILER("MatrixExponential::compute");
    }

    void computeExpTimesVector(const MatrixType& arg, const VectorType& v, VectorType& result, int vec_squarings)
    {
        START_PROFILER("MatrixExponential::computeExpTimesVector");
        START_PROFILER("MatrixExponential:computeExpTimesVector:computeUV");
        computeUV(arg, U, V, squarings); // Pade approximant is (U+V) / (-U+V)
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
        // number of squarings implemented via matrix-vector multiplications
        //    int vec_squarings = 8;
        // number of squarings implemented via matrix-matrix multiplications
        int mat_squarings = squarings - vec_squarings;

        for (int i = 0; i < mat_squarings; i++) {
            tmp2.noalias() = tmp * tmp;
            tmp = tmp2;
        }

        int two_pow_s = (int)std::pow(2, vec_squarings);
        v_tmp = v;
        for (int i = 0; i < two_pow_s; i++) {
            result.noalias() = tmp * v_tmp;
            v_tmp = result;
        }
        STOP_PROFILER("MatrixExponential:computeExpTimesVector:squaringVector");
        STOP_PROFILER("MatrixExponential::computeExpTimesVector");
    }

    void computeUV(const MatrixType& arg, MatrixType& U, MatrixType& V, int& squarings)
    {
        using std::frexp;
        using std::pow;
        const double l1norm = arg.cwiseAbs().colwise().sum().maxCoeff();
        squarings = 0;
        if (l1norm < 1.495585217958292e-002) {
            matrix_exp_pade3(arg, U, V);
        } else if (l1norm < 2.539398330063230e-001) {
            matrix_exp_pade5(arg, U, V);
        } else if (l1norm < 9.504178996162932e-001) {
            matrix_exp_pade7(arg, U, V);
        } else if (l1norm < 2.097847961257068e+000) {
            matrix_exp_pade9(arg, U, V);
        } else {
            const double maxnorm = 5.371920351148152;
            frexp(l1norm / maxnorm, &squarings);
            if (squarings < 0)
                squarings = 0;
            A_scaled = arg.unaryExpr(Eigen::internal::MatrixExponentialScalingOp<double>(squarings));
            START_PROFILER("MatrixExponential::matrix_exp_pade13");
            matrix_exp_pade13(A_scaled, U, V);
            STOP_PROFILER("MatrixExponential::matrix_exp_pade13");
        }
    }

    /** \brief Compute the (3,3)-Pad&eacute; approximant to the exponential.
   *
   *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
   *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
   */
    void matrix_exp_pade3(const MatrixType& A, MatrixType& U, MatrixType& V)
    {
        //    typedef typename Eigen::internal::NumTraits<typename traits<MatrixType>::Scalar>::Real RealScalar;
        typedef double RealScalar;
        const RealScalar b[] = { 120., 60., 12., 1. };
        A2.noalias() = A * A;
        tmp.noalias() = b[3] * A2 + b[1] * eye;
        U.noalias() = A * tmp;
        V.noalias() = b[2] * A2 + b[0] * eye;
    }

    /** \brief Compute the (5,5)-Pad&eacute; approximant to the exponential.
   *
   *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
   *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
   */
    void matrix_exp_pade5(const MatrixType& A, MatrixType& U, MatrixType& V)
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

    /** \brief Compute the (7,7)-Pad&eacute; approximant to the exponential.
   *
   *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
   *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
   */
    void matrix_exp_pade7(const MatrixType& A, MatrixType& U, MatrixType& V)
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

    /** \brief Compute the (9,9)-Pad&eacute; approximant to the exponential.
   *
   *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
   *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
   */
    void matrix_exp_pade9(const MatrixType& A, MatrixType& U, MatrixType& V)
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

    /** \brief Compute the (13,13)-Pad&eacute; approximant to the exponential.
   *
   *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
   *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
   */
    void matrix_exp_pade13(const MatrixType& A, MatrixType& U, MatrixType& V)
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
};

} // end namespace Eigen

#endif // EIGEN_MATRIX_EXPONENTIAL

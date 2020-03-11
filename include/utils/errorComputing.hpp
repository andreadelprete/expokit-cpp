#include <Eigen/Core>

/*
    Yes, this can be done better
    But it's a start since I needed this code in several different points
*/

// This function gets the result of the integration
// and computes the error in terms of Newton
Vector2d computeErrorIntegral24(Matrix<double, 24, 1>& x, Matrix<double, 24, 1>& xRef, double dt)
{
    const double K = 1e5;
    const double B = 1e2;
    VectorXd p0(12);
    p0 << -0.190000000000000002, 0.150050000000000017, 0.000053853008907090, -0.190000000000000002, -0.150050000000000017, 0.000053853008907090, 0.190000000000000002, 0.150050000000000017, 0.000053853008907090, 0.190000000000000002, -0.150050000000000017, 0.000053853008907090;

    VectorXd fInt = K * (p0 * dt - x.block(0, 0, 12, 1) - B * x.block(12, 0, 12, 1));
    VectorXd fIntRef = K * (p0 * dt - xRef.block(0, 0, 12, 1) - B * xRef.block(12, 0, 12, 1));

    VectorXd diff = fInt - fIntRef;

    Vector2d result;
    result(0) = diff.norm();
    result(1) = diff.lpNorm<Infinity>();

    return result;
}
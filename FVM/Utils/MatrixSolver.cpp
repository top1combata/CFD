#include "MatrixSolver.h"
#include <Eigen/Dense>


void relaxSystem(Matrix& A, Matrix& rhs, Matrix const& previousValue, Scalar relaxFactor)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());
    assert(rhs.rows() == previousValue.rows() && rhs.cols() == previousValue.cols());
    assert(relaxFactor > 0);
    
    for (Index idx = 0; idx < A.rows(); idx++)
    {
        rhs.row(idx) += (1/relaxFactor - 1) * A(idx, idx) * previousValue.row(idx);
        A(idx, idx) /= relaxFactor;
    }
}


Matrix solveSystem(Matrix& A, Matrix const& rhs)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());

    return A.partialPivLu().solve(rhs);
}
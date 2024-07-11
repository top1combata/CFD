#include "MatrixSolver.h"
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <cassert>

template<class AnyMatrix>
void relaxSystem(AnyMatrix& A, Matrix& rhs, Matrix const& previousValue, Scalar relaxFactor)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());
    assert(rhs.rows() == previousValue.rows() && rhs.cols() == previousValue.cols());
    assert(relaxFactor > 0);

    auto Adiag = A.diagonal();
    
    for (Index idx = 0; idx < A.rows(); idx++)
    {
        rhs.row(idx) += (1/relaxFactor - 1) * Adiag(idx) * previousValue.row(idx);
        Adiag(idx) /= relaxFactor;
    }
}
// explicit instantiation
template void relaxSystem<Matrix>(Matrix& A, Matrix& rhs, Matrix const& previousValue, Scalar relaxFactor);
template void relaxSystem<SparseMatrix>(SparseMatrix& A, Matrix& rhs, Matrix const& previousValue, Scalar relaxFactor);


Matrix solveSystem(Matrix& A, Matrix const& rhs)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());

    return A.partialPivLu().solve(rhs);
}


Matrix solveSystem(SparseMatrix& A, Matrix const& rhs)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());

    Eigen::BiCGSTAB<SparseMatrix> solver(A);
    return solver.solve(rhs);
}
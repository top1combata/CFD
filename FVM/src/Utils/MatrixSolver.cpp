#include "MatrixSolver.h"
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <cassert>


void relaxSystem(SparseMatrix& A, Matrix& rhs, Field<Scalar> const& previousValue, Scalar relaxFactor)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());
    assert(rhs.rows() == previousValue.rows() && rhs.cols() == 1);
    assert(relaxFactor > 0);

    auto Adiag = A.diagonal();
    
    for (Index idx = 0; idx < A.rows(); idx++)
    {
        rhs(idx) += (1/relaxFactor - 1) * Adiag(idx) * previousValue(idx);
        Adiag(idx) /= relaxFactor;
    }
}


void relaxSystem(SparseMatrix& A, Matrix& rhs, Field<Vector> const& previousValue, Scalar relaxFactor)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());
    assert(rhs.rows() == previousValue.rows() && rhs.cols() == 3);
    assert(relaxFactor > 0);

    auto Adiag = A.diagonal();
    
    for (Index idx = 0; idx < A.rows(); idx++)
    {
        rhs.row(idx) += (1/relaxFactor - 1) * Adiag(idx) * previousValue(idx).transpose();
        Adiag(idx) /= relaxFactor;
    }
}


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
    assert(solver.info() == Eigen::Success);

    Matrix result;
    result = solver.solve(rhs);

    assert(solver.info() == Eigen::Success);

    return result;
}


Matrix solveSystem(SparseMatrix& A, Matrix const& rhs, Matrix const& guess)
{
    assert(A.cols() == A.rows());
    assert(A.rows() == rhs.rows());
    assert(A.rows() == guess.rows());
    assert(rhs.cols() == guess.cols());

    Eigen::BiCGSTAB<SparseMatrix> solver(A);
    assert(solver.info() == Eigen::Success);

    Matrix result;
    result = solver.solveWithGuess(rhs, guess);

    assert(solver.info() == Eigen::Success);

    return result;
}
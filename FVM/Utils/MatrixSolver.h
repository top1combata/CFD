#pragma once

#include "Types.h"

// only for Matrix and SparseMatrix
template<class AnyMatrix>
void relaxSystem(AnyMatrix& A, Matrix& rhs, Matrix const& previousValue, Scalar relaxFactor);

Matrix solveSystem(Matrix& A, Matrix const& rhs);

Matrix solveSystem(SparseMatrix& A, Matrix const& rhs);
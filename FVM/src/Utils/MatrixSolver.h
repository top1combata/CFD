#pragma once

#include "Types.h"


// only for Matrix and SparseMatrix
void relaxSystem(SparseMatrix& A, Matrix& rhs, Field<Scalar> const& previousValue, Scalar relaxFactor);
void relaxSystem(SparseMatrix& A, Matrix& rhs, Field<Vector> const& previousValue, Scalar relaxFactor);

Matrix solveSystem(Matrix& A, Matrix const& rhs);

Matrix solveSystem(SparseMatrix& A, Matrix const& rhs);
Matrix solveSystem(SparseMatrix& A, Matrix const& rhs, Matrix const& guess);

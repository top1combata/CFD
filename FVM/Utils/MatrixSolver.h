#pragma once

#include "Types.h"

void relaxSystem(Matrix& A, Matrix& rhs, Matrix const& previousValue, Scalar relaxFactor = 1);

Matrix solveSystem(Matrix& A, Matrix const& rhs);
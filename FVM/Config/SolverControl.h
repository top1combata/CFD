#pragma once

#include "Utils/Types.h"

namespace Config
{
    constexpr Scalar uRelax = 0.3;
    constexpr Scalar pRelax = 0.1;

    constexpr Scalar uTolerance = 1e-3;
    constexpr Scalar pTolerance = 1e-3;
    constexpr Index maxIterations = 1000;
}
#pragma once

#include "Utils/Types.h"

namespace Config
{
    constexpr Scalar uRelax = 0.3;
    constexpr Scalar pRelax = 0.2;

    constexpr Scalar uTolerance = 1e-7;
    constexpr Scalar pTolerance = 1e-7;
    constexpr Index maxIterations = 10000;
}
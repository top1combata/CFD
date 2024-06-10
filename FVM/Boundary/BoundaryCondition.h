#pragma once

#include "Utils/Types.h"
#include <functional>

template<class T>
class BoundaryCondition
{
public:
    
    enum class BoundaryConditionType
    {
        FIXED_VALUE,
        FIXED_GRADIENT
    };

    T value;
    BoundaryConditionType type;
};

class Boundaries
{
    BoundaryCondition<Vector> uBoundary;
    BoundaryCondition<Scalar> pBoundary;
};


template<class T>
using BoundaryConditionGetter = std::function<BoundaryCondition<T>(Index)>;
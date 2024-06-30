#pragma once

#include "Utils/Types.h"
#include "Utils/TypesOperations.h"
#include <functional>


enum class BoundaryConditionType
{
    FIXED_VALUE,
    FIXED_GRADIENT
};

template<class T>
class BoundaryCondition
{
public:

    T value;
    BoundaryConditionType type;
};

class Boundaries
{
public:

    BoundaryCondition<Vector> uBoundary;
    BoundaryCondition<Scalar> pBoundary;

};


template<class T>
using BoundaryConditionGetter = std::function<BoundaryCondition<T>(Index)>;


template<class T>
BoundaryConditionGetter<T> zeroGrad()
{
    return [](Index) -> BoundaryCondition<T>
    {
        return {zero<T>(), BoundaryConditionType::FIXED_GRADIENT};
    };
}
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

    static BoundaryCondition fixedValue(T value)
    {
        return {value, BoundaryConditionType::FIXED_VALUE};
    }

    static BoundaryCondition fixedGradient(T value)
    {
        return {value, BoundaryConditionType::FIXED_GRADIENT};
    }
};

class Boundaries
{
public:

    BoundaryCondition<Vector> uBoundary;
    BoundaryCondition<Scalar> pBoundary;

    static Boundaries wall();
    static Boundaries movingWall(Vector velocity);
    static Boundaries outlet(Scalar pressure);
    static Boundaries inlet(Vector velocity);
};


template<class T>
using BoundaryConditionGetter = std::function<BoundaryCondition<T>(Index)>;


template<class T>
BoundaryConditionGetter<T> zeroGradGetter()
{
    return [](Index) -> BoundaryCondition<T>
    {
        static auto boundary = BoundaryCondition<T>::fixedGradient(zero<T>());
        return boundary;
    };
}

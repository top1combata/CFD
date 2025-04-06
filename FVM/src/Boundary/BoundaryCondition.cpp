#include "BoundaryCondition.h"


Boundaries Boundaries::wall()
{
    return 
    {
        .uBoundary{BoundaryCondition<Vector>::fixedValue({0, 0, 0})},
        .pBoundary{BoundaryCondition<Scalar>::fixedGradient(0)}
    };
}


Boundaries Boundaries::outlet(Scalar pressure)
{
    return 
    {
        .uBoundary{BoundaryCondition<Vector>::fixedGradient({0, 0, 0})},
        .pBoundary{BoundaryCondition<Scalar>::fixedValue(pressure)}
    };
}


Boundaries Boundaries::inlet(Vector velocity)
{
    return 
    {
        .uBoundary{BoundaryCondition<Vector>::fixedValue(velocity)},
        .pBoundary{BoundaryCondition<Scalar>::fixedGradient(0)}
    };
}
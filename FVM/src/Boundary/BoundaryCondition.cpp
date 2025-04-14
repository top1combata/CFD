#include "BoundaryCondition.h"


Boundaries Boundaries::wall()
{
    return 
    {
        .uBoundary{BoundaryCondition<Vector>::fixedValue(zero<Vector>())},
        .pBoundary{BoundaryCondition<Scalar>::fixedGradient(zero<Scalar>())}
    };
}


Boundaries Boundaries::movingWall(Vector velocity)
{
    return
    {
        .uBoundary{BoundaryCondition<Vector>::fixedValue(velocity)},
        .pBoundary{BoundaryCondition<Scalar>::fixedGradient(zero<Scalar>())}
    };
}


Boundaries Boundaries::outlet(Scalar pressure)
{
    return 
    {
        .uBoundary{BoundaryCondition<Vector>::fixedGradient(zero<Vector>())},
        .pBoundary{BoundaryCondition<Scalar>::fixedValue(pressure)}
    };
}


Boundaries Boundaries::inlet(Vector velocity)
{
    return 
    {
        .uBoundary{BoundaryCondition<Vector>::fixedValue(velocity)},
        .pBoundary{BoundaryCondition<Scalar>::fixedGradient(zero<Scalar>())}
    };
}

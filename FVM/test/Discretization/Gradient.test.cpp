#include "TestUtils.h"
#include "TestPrinting.h"

#include <Discretization/Interpolation.h>
#include <Mesh/2D/Structured/CartesianMesh2D.h>

#include <gtest/gtest.h>
#include <sstream>


TEST(TestGradient, GreenGauseScalarInteriorCell)
{
    constexpr Scalar dx = 0.1;
    constexpr Scalar dy = 0.3;

    CartesianMesh2D mesh(3, 3, 3 * dx, 3 * dy);

    Index centerCellIdx = 4;

    BoundaryConditionGetter<Vector> boundaryGetter = [](Index)
    {
        EXPECT_TRUE(false);
        return BoundaryCondition<Vector>{};
    };

    auto gradient = Interpolation::computeCellGradient
    (
        mesh, 
        centerCellIdx, 
        boundaryGetter, 
        Interpolation::Schemes::Gradient::GREEN_GAUSE
    );

    LinearCombination<Vector, Vector> expectedGradient;
    expectedGradient.terms = 
    {
        {Vector(0, -1 / (2 * dy), 0), 1},
        {Vector(0, 1 / (2 * dy), 0), 7},
        {Vector(-1 / (2 * dx), 0, 0), 3},
        {Vector(1 / (2 * dx), 0, 0), 5}
    };

    EXPECT_EQ(gradient, expectedGradient);
}

TEST(TestGradient, LeastSquareScalarInteriorCell)
{
    constexpr Scalar dx = 0.1;
    constexpr Scalar dy = 0.3;

    CartesianMesh2D mesh(3, 3, 3 * dx, 3 * dy);

    Index centerCellIdx = 4;

    BoundaryConditionGetter<Vector> boundaryGetter = [](Index)
    {
        EXPECT_TRUE(false);
        return BoundaryCondition<Vector>{};
    };

    auto gradient = Interpolation::computeCellGradient
    (
        mesh, 
        centerCellIdx, 
        boundaryGetter, 
        Interpolation::Schemes::Gradient::LEAST_SQAURE
    );

    LinearCombination<Vector, Vector> expectedGradient;
    expectedGradient.terms = 
    {
        {Vector(0, -1 / (2 * dy), 0), 1},
        {Vector(0, 1 / (2 * dy), 0), 7},
        {Vector(-1 / (2 * dx), 0, 0), 3},
        {Vector(1 / (2 * dx), 0, 0), 5}
    };

    EXPECT_EQ(gradient, expectedGradient);
}

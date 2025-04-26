#include <gtest/gtest.h>
#include "TestUtils.h"

#include "Solvers/SIMPLE/SimpleAlgorithm.h"
#include "Mesh/2D/Structured/CartesianMesh2D.h"
#include "Config/Config.h"


class PoiseuilleFixture : public testing::Test
{
protected:

    static constexpr Index nx = 10;
    static constexpr Index ny = 9;
    static constexpr Scalar lx = 1.0;
    static constexpr Scalar ly = 1.0;

    static constexpr Scalar inletPressure = 0.2;
    static constexpr Scalar outletPressure = 0.1;

    static constexpr Scalar density = 1e3;
    static constexpr Scalar viscosity = 1;

    static constexpr Scalar uRelax = 0.7;
    static constexpr Scalar pRelax = 0.3;
    
    static constexpr Scalar uTolerance = 1e-5;
    static constexpr Scalar pTolerance = 1e-7;
    static constexpr Index maxIterations = 10'000;

    static constexpr auto convectionScheme = Interpolation::Schemes::Convection::SOU;
    static constexpr auto gradientScheme = Interpolation::Schemes::Gradient::GREEN_GAUSE;


    PoiseuilleFixture() : m_mesh(nx, ny, lx, ly)
    {
        m_mesh.setLeftBoundary(Boundaries::outlet(inletPressure));
        m_mesh.setRightBoundary(Boundaries::outlet(outletPressure));
        m_mesh.setTopBoundary(Boundaries::wall());
        m_mesh.setBottomBoundary(Boundaries::wall());

        Config::density = density;
        Config::viscosity = viscosity;

        Config::uRelax = uRelax;
        Config::pRelax = pRelax;

        Config::uTolerance = uTolerance;
        Config::pTolerance = pTolerance;
        Config::maxIterations = maxIterations;

        Config::convectionScheme = convectionScheme;
        Config::gradientScheme = gradientScheme;
    }

    template<class Solver>
    void testSolver() const
    {
        Solver solver(m_mesh);
        solver.solve();

        ASSERT_TRUE(solver.isConverged());

        for (Index cellIdx = 0; cellIdx < m_mesh.getCellAmount(); cellIdx++)
        {
            Vector cellCenter = m_mesh.getCellCentroid(cellIdx);
            Scalar y = cellCenter.y() - ly / 2;
            Scalar x = cellCenter.x();
    
            Scalar expectedPressure = ((lx - x) * inletPressure + x * outletPressure) / lx;
            Scalar maxVelocity = (inletPressure - outletPressure) * ly * ly / (8 * viscosity * lx);
            Vector expectedVelocity = {maxVelocity * (1 - std::pow(y / (ly / 2), 2)), 0, 0};
    
            auto const& pressureField = solver.getPressure(0);
            auto const& velocityField = solver.getVelocity(0);

            Vector velocity = velocityField(cellIdx);
            Vector diff = expectedVelocity - velocity;
    
            EXPECT_NEAR(pressureField(cellIdx), expectedPressure, 1e-6);
            
            EXPECT_LE(diff.norm() / expectedVelocity.norm(), 0.1)
                << "value is " << to_string(velocity) << ", but expected " << to_string(expectedVelocity);
        }
    }

    CartesianMesh2D m_mesh;
};


TEST_F(PoiseuilleFixture, TestSimpleAlgorithm)
{
    testSolver<SimpleAlgorithm>();
}

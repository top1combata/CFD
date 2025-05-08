#include <Mesh/2D/Structured/CartesianMesh2D.h>
#include <Solvers/SIMPLE/SimpleAlgorithm.h>
#include <Config/Config.h>
#include <Visualization/PostProcessor.h>


int main(int argc, char** argv)
{
    Config::viscosity = 1;
    Config::uRelax = 0.3;
    Config::pRelax = 0.1;
    constexpr Scalar inletVelocity = 0.03;
    constexpr Scalar outletPressure = 0;

    CartesianMesh2D mesh(10, 10, 1.0, 0.2);
    mesh.setBottomBoundary(Boundaries::wall());
    mesh.setTopBoundary(Boundaries::wall());
    mesh.setRightBoundary(Boundaries::outlet(outletPressure));
    mesh.setLeftBoundary(Boundaries::inlet({inletVelocity, 0, 0}));

    SimpleAlgorithm solver(mesh);
    solver.solve();

    PostProcessor post(solver);
    post.show();

    return 0;
}

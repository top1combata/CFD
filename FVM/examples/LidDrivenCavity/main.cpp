#include <Mesh/2D/Structured/CartesianMesh2D.h>
#include <Solvers/SIMPLE/SimpleAlgorithm.h>
#include <Config/Config.h>
#include <Visualization/PostProcessor.h>


int main(int argc, char** argv)
{
    Config::viscosity = 1;
    Config::uRelax = 0.3;
    Config::pRelax = 0.1;
    Config::timeStep = 0.05;
    Config::timeBegin = 0;
    Config::timeEnd = 5;
    constexpr Scalar wallVelocity = 0.03;

    CartesianMesh2D mesh(50, 50, 1.0, 1.0);
    mesh.setBottomBoundary(Boundaries::outlet(0));
    mesh.setTopBoundary(Boundaries::movingWall({wallVelocity, 0, 0}));
    mesh.setRightBoundary(Boundaries::wall());
    mesh.setLeftBoundary(Boundaries::wall());

    SimpleAlgorithm SIMPLE(mesh);
    SIMPLE.setTransient(true);
    SIMPLE.solve();

    PostProcessor post(SIMPLE);
    post.show();

    return 0;
}

#include <Mesh/2D/Structured/CartesianMesh2D.h>
#include <Solvers/SIMPLE/SimpleAlgorithm.h>
#include <Config/Config.h>

#include <omp.h>
#include <fstream>

int main(int argc, char** argv)
{
    Config::viscosity = 1;
    Config::uRelax = 0.3;
    Config::pRelax = 0.1;
    constexpr Scalar wallVelocity = 0.03;

    CartesianMesh2D mesh(50, 50, 1.0, 1.0);
    mesh.setBottomBoundary(Boundaries::outlet(0));
    mesh.setTopBoundary(Boundaries::movingWall({wallVelocity, 0, 0}));
    mesh.setRightBoundary(Boundaries::wall());
    mesh.setLeftBoundary(Boundaries::wall());

    SimpleAlgorithm SIMPLE(mesh);
    SIMPLE.solve();

    std::ofstream U("U");
    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
        U << mesh.getCellCentroid(cellIdx).x() << ' ' << mesh.getCellCentroid(cellIdx).y() << ' ' << SIMPLE.getU()(cellIdx).x() << ' ' << SIMPLE.getU()(cellIdx).y() << '\n';

    return 0;
}

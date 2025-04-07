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
    constexpr Scalar inletVelocity = 0.03;
    constexpr Scalar outletPressure = 0;

    CartesianMesh2D mesh(15, 15, 1.0, 0.2);
    mesh.setBottomBoundary(Boundaries::wall());
    mesh.setTopBoundary(Boundaries::wall());
    mesh.setRightBoundary(Boundaries::outlet(outletPressure));
    mesh.setLeftBoundary(Boundaries::inlet({inletVelocity, 0, 0}));

    SimpleAlgorithm SIMPLE(mesh);
    SIMPLE.solve();

    std::ofstream U("U");
    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
        U << mesh.getCellCentroid(cellIdx).x() << ' ' << mesh.getCellCentroid(cellIdx).y() << ' ' << SIMPLE.getU()(cellIdx).x() << ' ' << SIMPLE.getU()(cellIdx).y() << '\n';

    return 0;
}

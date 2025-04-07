#include <Mesh/2D/Structured/CartesianMesh2D.h>
#include <Solvers/SIMPLE/SimpleAlgorithm.h>
#include <fstream>
#include <omp.h>
#include <Config/Config.h>


int main()
{
    Config::viscosity = 1;
    Config::uRelax = 0.3;
    Config::pRelax = 0.1;

    omp_set_num_threads(8);
    CartesianMesh2D mesh(15, 15, 1, 1);

    mesh.setBottomBoundary(Boundaries::wall());
    mesh.setTopBoundary(Boundaries::wall());
    mesh.setRightBoundary(Boundaries::outlet(0));
    mesh.setLeftBoundary(Boundaries::inlet({0.03, 0, 0}));

    SimpleAlgorithm SIMPLE(mesh);
    SIMPLE.solve();

    std::ofstream U("U");
    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
        U << mesh.getCellCentroid(cellIdx).x() << ' ' << mesh.getCellCentroid(cellIdx).y() << ' ' << SIMPLE.getU()(cellIdx).x() << ' ' << SIMPLE.getU()(cellIdx).y() << '\n';

    return 0;
}
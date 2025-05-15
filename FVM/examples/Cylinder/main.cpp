#include <Mesh/2D/Unstructured/PolyMesh2D.h>
#include <Solvers/SIMPLE/SimpleAlgorithm.h>
#include <Config/Config.h>
#include <Visualization/PostProcessor.h>

#include <sstream>


int main()
{
    Config::viscosity = 1;
    Config::uRelax = 0.3;
    Config::pRelax = 0.1;
    Config::timeStep = 0.005;
    Config::timeBegin = 0;
    Config::timeEnd = 1;
    Config::pTolerance = 1e-2;
    Config::maxIterations = 1000;
    Config::uSystemTolerance = 1e-3;
    Config::pSystemTolerance = 1e-3;

    PolyMesh2D mesh
    {
        std::stringstream
        (
            #include "cylinder_mesh.h"
        )
    };
    mesh.useNonOrthogonalCorrection = true;

    SimpleAlgorithm solver(mesh);
    solver.setTransient(true);
    solver.solve();

    PostProcessor post(solver);
    post.show();

    return 0;
}

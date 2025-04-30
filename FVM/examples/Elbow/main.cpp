#include <Mesh/2D/Unstructured/PolyMesh2D.h>
#include <Solvers/SIMPLE/SimpleAlgorithm.h>
#include <Config/Config.h>
#include <Visualization/PostProcessor.h>

#include <fstream>


int main(int argc, char** argv)
{
    if (argc < 2)
    {
        throw std::runtime_error("Need a mesh file");
    }

    Config::viscosity = 1;
    Config::uRelax = 0.3;
    Config::pRelax = 0.1;
    Config::timeStep = 0.0025;
    Config::timeBegin = 0;
    Config::timeEnd = 0.1;
    Config::pTolerance = 1e-4;
    Config::maxIterations = 1000;

    PolyMesh2D mesh{std::ifstream(argv[1])};
    mesh.useNonOrthogonalCorrection = true;

    SimpleAlgorithm SIMPLE(mesh);
    SIMPLE.setTransient(true);
    SIMPLE.solve();

    PostProcessor post(SIMPLE);
    post.show();

    return 0;
}

#include "Visualization/MeshDrawer.h"
#include "Mesh/2D/Unstructured/PolyMesh2D.h"

#include <fstream>
#include <filesystem>


int main(int argc, char** argv)
{
    if (argc < 2)
    {
        throw std::runtime_error("argument must be a mesh file");
    }

    std::string filename = argv[1];

    if (!std::filesystem::exists(filename))
    {
        throw std::runtime_error('\"' + filename + '\"' + ": no such file");
    }
    
    PolyMesh2D mesh{std::ifstream(filename)};
    MeshDrawer drawer(mesh);
    drawer.show();

    return 0;
}

#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"
#include "Boundary/BoundaryCondition.h"


class PolyMesh2D : public MeshBase
{
public:

    PolyMesh2D(std::istream& stream);
    PolyMesh2D(std::istream&& stream) : PolyMesh2D(stream) {}

    bool is2D() const override;
};

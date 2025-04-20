#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"
#include "Boundary/BoundaryCondition.h"

class CartesianMesh2D : public MeshBase
{
public:

    CartesianMesh2D(Index xSize, Index ySize, Scalar xLen, Scalar yLen);

    //boundaries
    void setTopBoundary(Boundaries boundaries);

    void setBottomBoundary(Boundaries boundaries);

    void setLeftBoundary(Boundaries boundaries);
    
    void setRightBoundary(Boundaries boundaries);

    bool is2D() const override;

protected:

    Index m_nx;
    Index m_ny;
};
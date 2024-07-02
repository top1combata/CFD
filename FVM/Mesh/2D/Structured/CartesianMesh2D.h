#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"
#include "Boundary/BoundaryCondition.h"

class CartesianMesh2D : public MeshBase
{
public:

    CartesianMesh2D(Index xSize, Index ySize, Scalar xLen, Scalar yLen);

    // cells
    Index getCellAmount() const override;

    List<Index> getCellNeighbours(Index cellIdx) const override;

    Vector getCellCentroid(Index cellIdx) const override;

    Scalar getCellVolume(Index cellIdx) const override;


    //faces
    Index getFaceAmount() const override;

    Vector getFaceCentroid(Index faceIdx) const override;

    bool isBoundaryFace(Index faceIdx) const override;

    Array<Index, 2> getFaceNeighbours(Index faceIdx) const override;

    Vector getFaceVector(Index faceIdx) const override;

    Boundaries getFaceBoundary(Index faceIdx) const override;

    List<CellFace> getCellFaces(Index cellIdx) const override;


    //boundaries
    void setTopBoundary(Boundaries boundaries);

    void setBottomBoundary(Boundaries boundaries);

    void setLeftBoundary(Boundaries boundaries);
    
    void setRightBoundary(Boundaries boundaries);

private:

    const Index m_x;
    const Index m_y;

    const Scalar m_xlen;
    const Scalar m_ylen;

    class CartesianBoundaries
    {
    public:
        Boundaries top, bottom, left, right;
    } 
    m_boundaries;

    Array<Index, 2> getLocalIndices(Index cellIdx) const;
    Index getGlobalIndex(Index xIdx, Index yIdx) const;
};
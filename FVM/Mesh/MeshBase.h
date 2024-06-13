#pragma once

#include "Utils/Types.h"
#include "Boundary/BoundaryCondition.h"


class MeshBase
{
public:

    // cells
    virtual Index getCellsAmount() const = 0;
    virtual List<Index> getNeighbourCells(Index cellIdx) const = 0;
    virtual Vector getCellCentroid(Index cellIdx) const = 0;
    virtual Scalar getCellVolume(Index cellIdx) const  = 0;

    //faces
    virtual Index getFacesAmount() const = 0;
    virtual Vector getFaceCentroid(Index faceIdx) const = 0;
    virtual bool isBoundaryFace(Index faceIdx) const = 0;
    virtual Array<Index, 2> getFaceNeighbours(Index faceIdx) const = 0;
    virtual Vector getFaceNormal(Index faceIdx, Index cellIdx) const = 0;
    virtual Boundaries getFaceBoundary(Index faceIdx) const = 0;

    class CellFace
    {
    public:
        Vector faceVector;
        Index faceIdx;
    };
    virtual List<CellFace> getCellFaces(Index cellIdx) const = 0;
};
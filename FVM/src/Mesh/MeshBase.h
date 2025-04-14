#pragma once

#include "Utils/Types.h"
#include "Boundary/BoundaryCondition.h"


class MeshBase
{
public:

    // cells
    Index  getCellAmount() const;

    Vector getCellCentroid(Index cellIdx) const;

    Scalar getCellVolume(Index cellIdx) const;

    //faces
    Index           getFaceAmount() const;

    Vector          getFaceCentroid(Index faceIdx) const;

    bool            isBoundaryFace(Index faceIdx) const;

    Array<Index, 2> getFaceNeighbors(Index faceIdx) const;

    Vector          getFaceVector(Index faceIdx) const;

    Boundaries      getFaceBoundary(Index faceIdx) const;

    Index           getFaceOwner(Index faceIdx) const;

    List<Index>     getCellFaces(Index cellIdx) const;


    bool useNonOrthogonalCorrection = false;

    virtual ~MeshBase() = default;

protected:

    List<Vector>         m_cellCentroids;
    List<Scalar>         m_cellVolumes;
    List<List<Index>>    m_cellFaces;
    List<Vector>         m_faceVectors;
    List<Vector>         m_faceCentroids;
    List<Array<Index,2>> m_faceNeighbors;

    HashMap<Index, Boundaries> m_boundariesMap;
};

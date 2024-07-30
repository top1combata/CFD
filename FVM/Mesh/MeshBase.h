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

    Array<Index, 2> getFaceNeighbours(Index faceIdx) const;

    Vector          getFaceVector(Index faceIdx) const;

    Boundaries      getFaceBoundary(Index faceIdx) const;

    Index           getFaceOwner(Index faceIdx) const;

    List<Index>     getCellFaces(Index cellIdx) const;

    virtual ~MeshBase() = default;

protected:

    List<Vector>         m_cell_centroids;
    List<Scalar>         m_cell_volumes;
    List<List<Index>>    m_cell_faces;
    List<Vector>         m_face_vectors;
    List<Vector>         m_face_centroids;
    List<Array<Index,2>> m_face_neighbors;

    HashMap<Index, Boundaries> m_boundaries_map;
};
#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"
#include "Boundary/BoundaryCondition.h"


class PolyMesh2D : public MeshBase
{
public:

    PolyMesh2D(std::istream& stream);
    
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


private:

    List<Vector> m_cell_centroids;

    List<Scalar> m_cell_volumes;

    List<Array<Vector,2>> m_faces;

    List<List<MeshBase::CellFace>> m_cell_faces;

    List<Array<Index,2>> m_face_neighbours;

    HashMap<Index, Boundaries> m_boundaries_map;
};
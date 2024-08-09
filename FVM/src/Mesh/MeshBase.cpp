#include "MeshBase.h"


Index MeshBase::getCellAmount() const
{
    return static_cast<Index>(m_cell_volumes.size());
}


Vector MeshBase::getCellCentroid(Index cellIdx) const
{
    return m_cell_centroids[cellIdx];
}


Scalar MeshBase::getCellVolume(Index cellIdx) const
{
    return m_cell_volumes[cellIdx];
}


Index MeshBase::getFaceAmount() const
{
    return static_cast<Index>(m_face_centroids.size());
}


Vector MeshBase::getFaceCentroid(Index faceIdx) const
{
    return m_face_centroids[faceIdx];
}


bool MeshBase::isBoundaryFace(Index faceIdx) const
{
    return -1 == m_face_neighbors[faceIdx][1];
}


Array<Index, 2> MeshBase::getFaceNeighbours(Index faceIdx) const
{
    return m_face_neighbors[faceIdx];
}


Vector MeshBase::getFaceVector(Index faceIdx) const
{
    return m_face_vectors[faceIdx];
}


Boundaries MeshBase::getFaceBoundary(Index faceIdx) const
{
    return m_boundaries_map.at(faceIdx);
}


Index MeshBase::getFaceOwner(Index faceIdx) const
{
    return m_face_neighbors[faceIdx][0];
}


List<Index> MeshBase::getCellFaces(Index cellIdx) const
{
    return m_cell_faces[cellIdx];
}
#include "MeshBase.h"


Index MeshBase::getCellAmount() const
{
    return static_cast<Index>(m_cellVolumes.size());
}


Vector MeshBase::getCellCentroid(Index cellIdx) const
{
    return m_cellCentroids[cellIdx];
}


Scalar MeshBase::getCellVolume(Index cellIdx) const
{
    return m_cellVolumes[cellIdx];
}


Index MeshBase::getFaceAmount() const
{
    return static_cast<Index>(m_faceCentroids.size());
}


Vector MeshBase::getFaceCentroid(Index faceIdx) const
{
    return m_faceCentroids[faceIdx];
}


bool MeshBase::isBoundaryFace(Index faceIdx) const
{
    return -1 == m_faceNeighbors[faceIdx][1];
}


Array<Index, 2> MeshBase::getFaceNeighbors(Index faceIdx) const
{
    return m_faceNeighbors[faceIdx];
}


Vector MeshBase::getFaceVector(Index faceIdx) const
{
    return m_faceVectors[faceIdx];
}


Boundaries MeshBase::getFaceBoundary(Index faceIdx) const
{
    return m_boundariesMap.at(faceIdx);
}


Index MeshBase::getFaceOwner(Index faceIdx) const
{
    return m_faceNeighbors[faceIdx][0];
}


List<Index> const& MeshBase::getCellFaces(Index cellIdx) const
{
    return m_cellFaces[cellIdx];
}

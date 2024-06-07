#include "CartesianMesh2D.h"

CartesianMesh2D::CartesianMesh2D(Index xSize, Index ySize, Scalar xLen, Scalar yLen) : 
    m_x(xSize), m_y(ySize), m_xlen(xLen), m_ylen(yLen) {}


Index CartesianMesh2D::getCellsAmount() const
{
    return m_x*m_y;
};


List<Index> CartesianMesh2D::getNeighbourCells(Index cellIdx) const
{
    auto [xIdx, yIdx] = getLocalIndices(cellIdx);
    List<Index> Neighbours;

    if (xIdx > 0)
        Neighbours.push_back(getGlobalIndex(xIdx-1, yIdx));
    if (xIdx < m_x-1)
        Neighbours.push_back(getGlobalIndex(xIdx+1, yIdx));
    if (yIdx > 0)
        Neighbours.push_back(getGlobalIndex(xIdx, yIdx-1));
    if (yIdx < m_y-1)
        Neighbours.push_back(getGlobalIndex(xIdx, yIdx+1));

    return Neighbours;
}


Vector CartesianMesh2D::getCellCentroid(Index cellIdx) const
{
    auto [xIdx, yIdx] = getLocalIndices(cellIdx);
    Scalar xCord = m_xlen/m_x*(2*xIdx+1)/2;
    Scalar yCord = m_ylen/m_y*(2*yIdx+1)/2;
    return {xCord, yCord, 0};
}


Index CartesianMesh2D::getFacesAmount() const
{
    return 2*m_x*m_y + m_x + m_y;
}


Vector CartesianMesh2D::getFaceCentroid(Index faceIdx) const
{
    Index idx  = faceIdx/4;
    Index side = faceIdx%4;
    Vector centroid = getCellCentroid(idx);

    Scalar dx = m_xlen/m_x, dy = m_ylen/m_y;
    switch (side)
    {
    case 0:
        centroid(1) -= dy/2;
        break;
    case 1:
        centroid(0) += dx/2;
        break;
    case 2:
        centroid(1) += dy/2;
        break;
    case 3:
        centroid(0) -= dx/2;
        break;
    }
    return centroid;
}


bool CartesianMesh2D::isBoundaryFace(Index faceIdx) const
{
    Index idx = faceIdx/4, side = faceIdx%4;
    auto [xIdx, yIdx] = getLocalIndices(idx);

    return xIdx == 0 && side == 3 || xIdx == m_x-1 && side == 1 || yIdx == 0 && side == 0 || yIdx == m_y-1 && side == 2;
}


Boundaries CartesianMesh2D::getFaceBoundary(Index faceIdx) const
{
    Index side = faceIdx%4;
    switch (side)
    {
    case 0:
        return m_boundaries.bottom;
    case 1:
        return m_boundaries.right;
    case 2:
        return m_boundaries.top;
    case 3:
        return m_boundaries.left;
    }
}


List<MeshBase::CellFace> CartesianMesh2D::getCellFaces(Index cellIdx) const
{
    Scalar dx = m_xlen/m_x, dy = m_ylen/m_y;
    return {{{0,-dy,0}, 4*cellIdx}, {{dx,0,0}, 4*cellIdx+1}, {{0,dy,0}, 4*cellIdx+2}, {{-dx,0,0}, 4*cellIdx+3}};
}


Array<Index, 2> CartesianMesh2D::getLocalIndices(Index cellIdx) const
{
    return {cellIdx%m_x, cellIdx/m_x};
}


Index CartesianMesh2D::getGlobalIndex(Index xIdx, Index yIdx) const
{
    return yIdx*m_x + xIdx;
}
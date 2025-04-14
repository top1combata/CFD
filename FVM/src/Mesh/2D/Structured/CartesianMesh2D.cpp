#include "CartesianMesh2D.h"


CartesianMesh2D::CartesianMesh2D(Index xSize, Index ySize, Scalar xLen, Scalar yLen) :
    m_nx(xSize),
    m_ny(ySize)
{
    Scalar dx = xLen / xSize;
    Scalar dy = yLen / ySize;

    Index totalCells = xSize*ySize;
    Index totalFaces = 2*xSize*ySize + xSize + ySize;
    
    m_cellVolumes = List<Scalar>(totalCells, dx*dy);
    m_cellCentroids.resize(totalCells);
    m_cellFaces.resize(totalCells);
    m_faceNeighbors = List<Array<Index, 2>>(totalFaces, {-1,-1});
    m_faceCentroids.resize(totalFaces);
    m_faceVectors.resize(totalFaces);

    for (Index y = 0; y < ySize; y++)
    {
        for (Index x = 0; x < xSize; x++)
        {
            Index cellIdx = y*xSize + x;
            
            Index face1 = y*(2*xSize+1) + x;
            Index face2 = face1 + xSize;
            Index face3 = face1 + 2*xSize + 1;
            Index face4 = face2 + 1;

            m_cellFaces[cellIdx] = {face1, face2, face3, face4};

            Vector cellCentroid = {dx*(x+0.5), dy*(y+0.5), 0};
            m_cellCentroids[cellIdx] = cellCentroid;
            m_faceCentroids[face1] = cellCentroid + Vector{0,-dy/2,0};
            m_faceCentroids[face2] = cellCentroid + Vector{-dx/2,0,0};
            m_faceCentroids[face3] = cellCentroid + Vector{0,dy/2,0};
            m_faceCentroids[face4] = cellCentroid + Vector{dx/2,0,0};

            for (Index faceIdx : m_cellFaces[cellIdx])
            {
                if (m_faceNeighbors[faceIdx][0] == -1)
                {
                    m_faceNeighbors[faceIdx][0] = cellIdx;
                }
                else
                {
                    m_faceNeighbors[faceIdx][1] = cellIdx;
                }

                if (faceIdx % (2*xSize+1) < xSize)
                {
                    m_faceVectors[faceIdx] = {0,dx,0};
                }
                else
                {
                    m_faceVectors[faceIdx] = {dy,0,0};
                }
                
                if (faceIdx < xSize || (faceIdx-xSize) % (2*xSize+1) == 0)
                {
                    m_faceVectors[faceIdx] *= -1;
                }
            }
        }
    }
}


void CartesianMesh2D::setTopBoundary(Boundaries boundaries)
{
    for (Index faceIdx = 0; faceIdx < m_nx; faceIdx++)
    {
        m_boundariesMap[faceIdx] = boundaries;
    }
}


void CartesianMesh2D::setBottomBoundary(Boundaries boundaries)
{
    Index totalFaces = getFaceAmount();
    for (Index faceIdx = totalFaces-1; faceIdx >= totalFaces - m_nx; faceIdx--)
    {
        m_boundariesMap[faceIdx] = boundaries;
    }
}


void CartesianMesh2D::setLeftBoundary(Boundaries boundaries)
{
    for (Index faceIdx = m_nx; faceIdx < getFaceAmount(); faceIdx += 2*m_nx+1)
    {
        m_boundariesMap[faceIdx] = boundaries;
    }
}   


void CartesianMesh2D::setRightBoundary(Boundaries boundaries)
{
    for (Index faceIdx = 2*m_nx; faceIdx < getFaceAmount(); faceIdx += 2*m_nx+1)
    {
        m_boundariesMap[faceIdx] = boundaries;
    }
}

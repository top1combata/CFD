#include "CartesianMesh2D.h"

CartesianMesh2D::CartesianMesh2D(Index xSize, Index ySize, Scalar xLen, Scalar yLen) :
    m_nx(xSize),
    m_ny(ySize)
{
    Scalar dx = xLen / xSize;
    Scalar dy = yLen / ySize;

    Index totalCells = xSize*ySize;
    Index totalFaces = 2*xSize*ySize + xSize + ySize;
    
    m_cell_volumes = List<Scalar>(totalCells, dx*dy);
    m_cell_centroids.resize(totalCells);
    m_cell_faces.resize(totalCells);
    m_face_neighbors = List<Array<Index, 2>>(totalFaces, {-1,-1});
    m_face_centroids.resize(totalFaces);
    m_face_vectors.resize(totalFaces);

    for (Index y = 0; y < ySize; y++)
    {
        for (Index x = 0; x < xSize; x++)
        {
            Index cellIdx = y*xSize + x;
            
            Index face1 = y*(2*xSize+1) + x;
            Index face2 = face1 + xSize;
            Index face3 = face1 + 2*xSize + 1;
            Index face4 = face2 + 1;

            m_cell_faces[cellIdx] = {face1, face2, face3, face4};

            Vector cellCentroid = {dx*(x+0.5), dy*(y+0.5), 0};
            m_cell_centroids[cellIdx] = cellCentroid;
            m_face_centroids[face1] = cellCentroid + Vector{0,-dy/2,0};
            m_face_centroids[face2] = cellCentroid + Vector{-dx/2,0,0};
            m_face_centroids[face3] = cellCentroid + Vector{0,dy/2,0};
            m_face_centroids[face4] = cellCentroid + Vector{dx/2,0,0};

            for (Index faceIdx : m_cell_faces[cellIdx])
            {
                if (m_face_neighbors[faceIdx][0] == -1)
                    m_face_neighbors[faceIdx][0] = cellIdx;
                else
                    m_face_neighbors[faceIdx][1] = cellIdx;

                if (faceIdx % (2*xSize+1) < xSize)
                    m_face_vectors[faceIdx] = {0,dx,0};
                else
                    m_face_vectors[faceIdx] = {dy,0,0};
                
                if (faceIdx < xSize || (faceIdx-xSize) % (2*xSize+1) == 0)
                    m_face_vectors[faceIdx] *= -1;
            }
        }
    }
}


void CartesianMesh2D::setTopBoundary(Boundaries boundaries)
{
    for (Index faceIdx = 0; faceIdx < m_nx; faceIdx++)
        m_boundaries_map[faceIdx] = boundaries;
}


void CartesianMesh2D::setBottomBoundary(Boundaries boundaries)
{
    Index totalFaces = getFaceAmount();
    for (Index faceIdx = totalFaces-1; faceIdx >= totalFaces - m_nx; faceIdx--)
        m_boundaries_map[faceIdx] = boundaries;
}


void CartesianMesh2D::setLeftBoundary(Boundaries boundaries)
{
    for (Index faceIdx = m_nx; faceIdx < getFaceAmount(); faceIdx += 2*m_nx+1)
        m_boundaries_map[faceIdx] = boundaries;
}   


void CartesianMesh2D::setRightBoundary(Boundaries boundaries)
{
    for (Index faceIdx = 2*m_nx; faceIdx < getFaceAmount(); faceIdx += 2*m_nx+1)
        m_boundaries_map[faceIdx] = boundaries;
}
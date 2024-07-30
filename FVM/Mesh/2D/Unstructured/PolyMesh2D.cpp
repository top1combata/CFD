#include "PolyMesh2D.h"
#include "Parse.h"


PolyMesh2D::PolyMesh2D(std::istream& stream)
{
    List<Vector> vertices;
    List<Array<Index,2>> faces;
    
    parse(stream, vertices, faces, m_cell_faces, m_boundaries_map);

    Index totalCells = m_cell_faces.size();
    Index totalFaces = faces.size();

    m_cell_centroids.resize(totalCells);
    m_cell_volumes.resize(totalCells);
    m_face_vectors.resize(totalFaces);
    m_face_centroids.resize(totalFaces);
    m_face_neighbors.resize(totalFaces, {-1,-1});

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        Vector sum = {0,0,0};
        for (Index faceIdx : m_cell_faces[cellIdx])
        {
            // assigning the neighbour cell
            if (-1 == m_face_neighbors[faceIdx][0])
                m_face_neighbors[faceIdx][0] = cellIdx;
            else
            {
                m_face_neighbors[faceIdx][1] = cellIdx;
                if (cellIdx < m_face_neighbors[faceIdx][0])
                    std::swap(m_face_neighbors[faceIdx][0], m_face_neighbors[faceIdx][1]);  
            }

            Index idx1 = faces[faceIdx][0];
            Index idx2 = faces[faceIdx][1];
            Vector v1 = vertices.at(idx1);
            Vector v2 = vertices.at(idx2);

            m_face_centroids[faceIdx] = (v1+v2)/2;
            m_face_vectors[faceIdx] = {(v2-v1)[1], -(v2-v1)[0], 0};
            sum += v1+v2;
        }
        m_cell_centroids.at(cellIdx) = sum / (2*m_cell_faces[cellIdx].size());

        auto triangleArea = [](Vector v1, Vector v2, Vector v3) -> Scalar
        {
            Vector r1 = v2-v1;
            Vector r2 = v3-v1;
            return abs(r1(0)*r2(1) - r1(1)*(r2(0)))/2;
        };

        m_cell_volumes[cellIdx] = 0;
        for (Index faceIdx : m_cell_faces[cellIdx])
        {
            Vector v1 = vertices.at(faces[faceIdx][0]);
            Vector v2 = vertices.at(faces[faceIdx][1]);
            Vector v3 = m_cell_centroids[cellIdx];

            m_cell_volumes[cellIdx] += triangleArea(v1,v2,v3);
        }
    }
    // correcting face vectors to neighbors
    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Vector vertex = vertices.at(faces[faceIdx][0]);
        Index owner = getFaceOwner(faceIdx);
        Vector r = m_cell_centroids.at(owner) - vertex;
        
        if (getFaceVector(faceIdx).dot(r) > 0)
            m_face_vectors[faceIdx] *= -1;
    }
}
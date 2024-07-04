#include "PolyMesh2D.h"
#include "Parse.h"

Index PolyMesh2D::getCellAmount() const
{
    return static_cast<Index>(m_cell_centroids.size());
}


List<Index> PolyMesh2D::getCellNeighbours(Index cellIdx) const
{
    List<Index> neighbours;
    neighbours.reserve(m_cell_faces[cellIdx].size());

    for (auto [_, faceIdx] : m_cell_faces[cellIdx])
    {
        auto [cell1, cell2] = getFaceNeighbours(faceIdx);
        if (isBoundaryFace(faceIdx) == false)
            neighbours.push_back(cell1 == cellIdx ? cell2 : cell1);
    }
    
    return neighbours;
}


Vector PolyMesh2D::getCellCentroid(Index cellIdx) const
{
    return m_cell_centroids[cellIdx];
}


Scalar PolyMesh2D::getCellVolume(Index cellIdx) const
{
    return m_cell_volumes[cellIdx];
}


Index PolyMesh2D::getFaceAmount() const
{
    return static_cast<Index>(m_faces.size());
}


Vector PolyMesh2D::getFaceCentroid(Index faceIdx) const
{
    return (m_faces[faceIdx][0] + m_faces[faceIdx][1]) / 2;
}


bool PolyMesh2D::isBoundaryFace(Index faceIdx) const
{
    return m_boundaries_map.contains(faceIdx);
}


Array<Index, 2> PolyMesh2D::getFaceNeighbours(Index faceIdx) const
{
    return m_face_neighbours[faceIdx];
}


Vector PolyMesh2D::getFaceVector(Index faceIdx) const
{
    Vector vec = m_faces[faceIdx][1] - m_faces[faceIdx][0];
    return {-vec[1], vec[0], 0};
    
}


Boundaries PolyMesh2D::getFaceBoundary(Index faceIdx) const
{
    return m_boundaries_map.at(faceIdx);
}


List<MeshBase::CellFace> PolyMesh2D::getCellFaces(Index cellIdx) const
{
    return m_cell_faces[cellIdx];
}



PolyMesh2D::PolyMesh2D(std::istream& stream)
{
    List<Vector> vertices;
    List<Array<Index,2>> faces;
    List<List<Index>> cells;
    
    parse(stream, vertices, faces, cells, m_boundaries_map);

    Index totalCells = cells.size();
    Index totalFaces = faces.size();
    m_cell_centroids.resize(totalCells);
    m_cell_volumes.resize(totalCells);
    m_cell_faces.resize(totalCells);
    m_face_neighbours.resize(totalFaces);
    m_faces.resize(totalFaces);

    // sorting indices
    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
        if (faces.at(faceIdx)[0] > faces.at(faceIdx)[1])
            std::swap(faces.at(faceIdx)[0], faces.at(faceIdx)[1]);

    auto pairHash = [](std::pair<Index, Index> pair)
    {
        std::hash<Index> h;
        return h(pair.first) ^ h(pair.second);
    };
    HashMap<std::pair<Index, Index>, List<Index>, decltype(pairHash)> faceCellMap;

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        Vector sum = {0,0,0};
        for (Index faceIdx : cells.at(cellIdx))
        {
            Index idx1 = faces.at(faceIdx).at(0);
            Index idx2 = faces.at(faceIdx).at(1);
            Vector v1 = vertices.at(idx1);
            Vector v2 = vertices.at(idx2);
            
            m_faces.at(faceIdx) = {v1,v2};
            faceCellMap[{idx1, idx2}].push_back(cellIdx);
            sum += v1+v2;
        }
        m_cell_centroids.at(cellIdx) = sum / (2*cells.at(cellIdx).size());

        auto triangleArea = [](Vector v1, Vector v2, Vector v3) -> Scalar
        {
            Vector r1 = v2-v1;
            Vector r2 = v3-v1;
            return abs(r1(0)*r2(1) - r1(1)*(r2(0)))/2;
        };

        m_cell_volumes[cellIdx] = 0;
        for (Index faceIdx : cells.at(cellIdx))
        {
            if (getFaceVector(faceIdx).dot(m_faces[faceIdx][0] - m_cell_centroids[cellIdx]) < 0)
                std::swap(m_faces[faceIdx][0], m_faces[faceIdx][1]);

            m_cell_faces.at(cellIdx).push_back({getFaceVector(faceIdx), faceIdx});
            m_cell_volumes[cellIdx] += triangleArea(m_cell_centroids[cellIdx], m_faces[faceIdx][0], m_faces[faceIdx][1]);
        }
    }

    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        List<Index>& cellIds = faceCellMap[{faces[faceIdx][0], faces[faceIdx][1]}];
        if (cellIds.size() == 2)
            m_face_neighbours[faceIdx] = {cellIds[0], cellIds[1]};
        else
            m_face_neighbours[faceIdx] = {cellIds[0], -1};
    }
}
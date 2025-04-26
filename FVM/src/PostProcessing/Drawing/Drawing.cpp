#include "Drawing.h"


sf::View getInitializeView(MeshBase const& mesh)
{
    auto limits = std::numeric_limits<Scalar>();
    Scalar minX = limits.max();
    Scalar minY = limits.max();
    Scalar maxX = limits.min();
    Scalar maxY = limits.min();

    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
    {
        Vector cellCentroid = mesh.getCellCentroid(cellIdx);

        minX = std::min(minX, cellCentroid.x());
        maxX = std::max(maxX, cellCentroid.x());
        minY = std::min(minY, cellCentroid.y());
        maxY = std::max(maxY, cellCentroid.y());
    }

    sf::View view;
    view.setCenter(sf::Vector2f((minX + maxX) / 2, (minY + maxY) / 2));
    float size = std::max(maxX - minX, maxY - minY);
    view.setSize(1.25f * sf::Vector2f(size, size));
    return view;
}


List<sf::ConvexShape> getMeshCells(MeshBase const& mesh)
{    
    auto getCellPolygon = [&mesh](Index cellIdx)
    {
        sf::ConvexShape polygon(mesh.getCellFaces(cellIdx).size());

        Index pointIdx = 0;
        for (Index faceIdx : mesh.getCellFaces(cellIdx))
        {
            Vector faceVector = mesh.getFaceVector(faceIdx);
            if (mesh.getFaceOwner(faceIdx) != cellIdx)
            {
                faceVector *= -1;
            }

            Vector cellVertex = 
                mesh.getFaceCentroid(faceIdx) - 
                Scalar(0.5) * Vector{faceVector.y(), -faceVector.x(), 0};

            polygon.setPoint(pointIdx++, sf::Vector2f(cellVertex.x(), cellVertex.y()));
        }
        return polygon;
    };

    List<sf::ConvexShape> cells(mesh.getCellAmount());
    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
    {
        cells[cellIdx] = getCellPolygon(cellIdx);
    }
    return cells;
}


List<Line> getMeshFaces(MeshBase const& mesh)
{
    auto getFace = [&mesh](Index faceIdx)
    {
        Vector faceVector = mesh.getFaceVector(faceIdx);
        Vector orthogonal{faceVector.y(), -faceVector.x(), 0};
        Vector begin = mesh.getFaceCentroid(faceIdx) - orthogonal / 2;
        Vector end = begin + orthogonal;
        
        return Line(begin, end);
    };

    List<Line> faces;
    faces.reserve(mesh.getFaceAmount());
    for (Index faceIdx = 0; faceIdx < mesh.getFaceAmount(); faceIdx++)
    {
        faces.emplace_back(getFace(faceIdx));
    }
    return faces;
}


List<Arrow> getVelocityField(SolverBase const& solver)
{
    auto& field = solver.getVelocity(0);
    auto const& mesh = solver.getMesh();
    
    List<Arrow> arrows;
    arrows.reserve(mesh.getCellAmount());
    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
    {
        arrows.emplace_back
        (
            mesh.getCellCentroid(cellIdx),
            field(cellIdx)
        );
    }

    return arrows;
}

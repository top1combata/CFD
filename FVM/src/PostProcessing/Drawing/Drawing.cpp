#include "Drawing.h"

#include "Mesh/Geometry.h"


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
    constexpr sf::Color FACE_COLOR = {150, 150, 150};

    auto getFace = [&](Index faceIdx)
    {
        Vector faceVector = mesh.getFaceVector(faceIdx);
        Vector orthogonal{faceVector.y(), -faceVector.x(), 0};
        Vector begin = mesh.getFaceCentroid(faceIdx) - orthogonal / 2;
        Vector end = begin + orthogonal;
        
        Line line(begin, end);
        line.setColor(FACE_COLOR);
        return line;
    };

    List<Line> faces;
    faces.reserve(mesh.getFaceAmount());
    for (Index faceIdx = 0; faceIdx < mesh.getFaceAmount(); faceIdx++)
    {
        faces.emplace_back(getFace(faceIdx));
    }
    return faces;
}


Vector getCellBoundedVector(MeshBase const& mesh, Index cellIdx, Vector direction)
{
    direction.normalize();
    Scalar distance = std::numeric_limits<Scalar>::max();

    for (Index faceIdx : mesh.getCellFaces(cellIdx))
    {
        Scalar distanceToFace = Geometry::distanceCellToFaceInDirection(mesh, cellIdx, faceIdx, direction);
        distance = std::min(distance, std::abs(distanceToFace));
    }

    return distance * direction;
}


List<Arrow> getVelocityField(SolverBase const& solver, Index timePointIdx, bool useAbsoluteVectorLength)
{
    auto& field = solver.getVelocity(timePointIdx);
    auto const& mesh = solver.getMesh();
    
    List<Arrow> arrows;
    arrows.reserve(mesh.getCellAmount());
    for (Index cellIdx = 0; cellIdx < mesh.getCellAmount(); cellIdx++)
    {
        Vector direction = field(cellIdx);
        if (!useAbsoluteVectorLength)
        {
            direction = getCellBoundedVector(mesh, cellIdx, direction);
        }

        arrows.emplace_back
        (
            mesh.getCellCentroid(cellIdx),
            direction
        );
    }

    return arrows;
}


sf::Color getGradientColor(Scalar t)
{
    t = 1 - t;
    if (t < Scalar(1) / 3)
    {
        return sf::Color(255, 3 * t * 255, 0);
    }
    if (t < Scalar(2) / 3)
    {
        return sf::Color(255 - (3 * t - 1) * 255, 255, 0);
    }
    return sf::Color(0, 255, (3 * t - 2) * 255);
}


void heatMapColor(List<sf::ConvexShape>& cells, std::function<Scalar(Index)> const& fieldGetter)
{
    Index totalCells = cells.size();
    
    Scalar minValue = std::numeric_limits<Scalar>::max();
    Scalar maxValue = std::numeric_limits<Scalar>::min();
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        Scalar fieldValue = fieldGetter(cellIdx);
        minValue = std::min(minValue, fieldValue);
        maxValue = std::max(maxValue, fieldValue);
    }

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        Scalar fieldValue = fieldGetter(cellIdx);
        Scalar t = (fieldValue - minValue) / (maxValue - minValue);
        cells[cellIdx].setFillColor(getGradientColor(t));
    }
}

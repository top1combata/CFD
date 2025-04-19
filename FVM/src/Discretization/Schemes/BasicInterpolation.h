#pragma once

#include "Discretization/LinearCombination.h"
#include "Mesh/Geometry.h"


namespace Interpolation
{

template<class T>
LinearCombination<T, Scalar> valueOnFace
(
    MeshBase const& mesh, 
    Index faceIdx, 
    BoundaryConditionGetter<T> const& boundaries
)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        auto [boundaryValue, boundaryType] = boundaries(faceIdx);

        using enum BoundaryConditionType;

        switch (boundaryType)
        {
        case FIXED_VALUE:
            return LinearCombination<T, Scalar>(boundaryValue);

        case FIXED_GRADIENT:
            Index cellIdx = mesh.getFaceNeighbors(faceIdx).front();
            Scalar dist = Geometry::distanceCellToFace(mesh, cellIdx, faceIdx);

            LinearCombination<T> result = {{1, cellIdx}};
            result += dist*boundaryValue;
            return  result;
        }

    }

    // general case
    auto [cellFromIdx, cellToIdx] = mesh.getFaceNeighbors(faceIdx);

    Scalar coeffFrom, coeffTo;

    if (mesh.useNonOrthogonalCorrection)
    {
        Scalar distanceBetweenCells = Geometry::distanceCellToCell(mesh, cellFromIdx, cellToIdx);
        Vector unitDirection = Geometry::cellToCellUnitVector(mesh, cellFromIdx, cellToIdx);
        Scalar distanceToFace = Geometry::distanceCellToFaceInDirection(mesh, cellFromIdx, faceIdx, unitDirection);
    
        coeffFrom = distanceToFace / distanceBetweenCells;
        coeffTo = 1 - coeffFrom;
    }
    else 
    {
        Scalar distFrom = Geometry::distanceCellToFace(mesh, cellFromIdx, faceIdx);
        Scalar distTo = Geometry::distanceCellToFace(mesh, cellToIdx, faceIdx);
    
        coeffFrom = 1 - distFrom/ (distFrom + distTo);
        coeffTo = 1 - coeffFrom;
    }
    return {{coeffFrom, cellFromIdx}, {coeffTo, cellToIdx}};
}



} // namespace Interpolation

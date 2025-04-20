#pragma once

#include "InterpolationSchemes.h"
#include "GradientSchemes.h"


namespace Interpolation
{

template<class T>
LinearCombination<T, Vector> computeCellGradient
(
    MeshBase const& mesh, 
    Index cellIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Schemes::Gradient::Type schemeType
)
{
    switch (schemeType)
    {
        case Schemes::Gradient::GREEN_GAUSE:
            return Schemes::Gradient::greenGauseGradientImpl
            (
                mesh, cellIdx, boundaries
            );

        case Schemes::Gradient::LEAST_SQAURE:
            return Schemes::Gradient::leastSquareGradientImpl
            (
                mesh, cellIdx, boundaries
            );
    }
}


// Warning ! Not for boundary faces
template<class T>
LinearCombination<T, Vector> computeAverageFaceGradient
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Schemes::Gradient::Type schemeType
)
{
    auto [ownerIdx, neighborIdx] = mesh.getFaceNeighbors(faceIdx);

    auto ownerCellGradient = computeCellGradient(mesh, ownerIdx, boundaries, schemeType);
    auto neighborCellGradient = computeCellGradient(mesh, neighborIdx, boundaries, schemeType);

    Vector faceCentroid = mesh.getFaceCentroid(faceIdx);
    Scalar ownerDistance = (faceCentroid - mesh.getCellCentroid(ownerIdx)).norm();
    Scalar neighborDistance = (faceCentroid - mesh.getCellCentroid(neighborIdx)).norm();

    Scalar ownerWeight = ownerDistance / (ownerDistance + neighborDistance);
    Scalar neighborWeight = 1 - ownerWeight;
    
    return ownerWeight * ownerCellGradient + neighborWeight * neighborCellGradient;
}


// Warning ! Not for boundary faces
template<class T>
LinearCombination<T, Vector> computeFaceGradient
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Schemes::Gradient::Type schemeType
)
{
    auto avgFaceGradient = computeAverageFaceGradient(mesh, faceIdx, boundaries, schemeType);
    
    if (!mesh.useNonOrthogonalCorrection)
    {
        return avgFaceGradient;
    }

    auto [ownerIdx, neighborIdx] = mesh.getFaceNeighbors(faceIdx);
    Vector unitDirection = Geometry::cellToCellUnitVector(mesh, ownerIdx, neighborIdx);
    Scalar distanceBetweenCells = Geometry::distanceCellToCell(mesh, ownerIdx, neighborIdx);

    LinearCombination<T, Scalar> finiteDifference = 
    {
        { 1 / distanceBetweenCells, ownerIdx},
        {-1 / distanceBetweenCells, neighborIdx}
    };

    auto nonOrthogonalCorrection = (finiteDifference - avgFaceGradient.dot(unitDirection)) * unitDirection;

    return avgFaceGradient + nonOrthogonalCorrection;
}


template<class T>
LinearCombination<T, Scalar> computeFaceNormalGradient
(
    MeshBase const& mesh,
    Index cellFromIdx, 
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
        case FIXED_GRADIENT:
            return LinearCombination<T, Scalar>(boundaryValue);
        
        case FIXED_VALUE:
            Scalar dist = Geometry::distanceCellToFace(mesh, cellFromIdx, faceIdx);

            LinearCombination<T> result = {{-1/dist, cellFromIdx}};
            result += boundaryValue / dist;
            return result;
        }
    }

    // General case
    Index cellToIdx;
    {
        auto [tmpIdx1, tmpIdx2] = mesh.getFaceNeighbors(faceIdx);
        cellToIdx = tmpIdx1 + tmpIdx2 - cellFromIdx;
    }

    Scalar distanceBetweenCells = Geometry::distanceCellToCell(mesh, cellFromIdx, cellToIdx);
    LinearCombination<T, Scalar> finiteDifference = 
    {
        {1 / distanceBetweenCells, cellToIdx},
        {-1 / distanceBetweenCells, cellFromIdx},
    };

    if (!mesh.useNonOrthogonalCorrection)
    {
        return finiteDifference;
    }

    // Least square cause it better for non-orthogonal meshes
    auto gradientScheme = Schemes::Gradient::LEAST_SQAURE;
    auto avgFaceGradient = computeAverageFaceGradient(mesh, faceIdx, boundaries, gradientScheme);

    Vector faceVector = mesh.getFaceVector(faceIdx).normalized();
    if (mesh.getFaceOwner(faceIdx) != cellFromIdx)
    {
        faceVector *= -1;
    }

    // Using over-relaxed approach
    Vector unitDirection = Geometry::cellToCellUnitVector(mesh, cellFromIdx, cellToIdx);
    Vector orthogonalComponent = faceVector.norm() / faceVector.dot(unitDirection) * unitDirection;
    Vector nonOrthogonalComponent = faceVector - orthogonalComponent;

    return 
    (
        finiteDifference * orthogonalComponent.norm() +
        avgFaceGradient.dot(nonOrthogonalComponent)
    );
}

} // namespace Interpolation::Schemes::Gradient

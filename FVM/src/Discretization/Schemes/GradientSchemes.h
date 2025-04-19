#pragma once

#include "Mesh/MeshBase.h"
#include "Mesh/Geometry.h"
#include "Discretization/LinearCombination.h"
#include "BasicInterpolation.h"

#include <Eigen/LU>


namespace Interpolation::Schemes::Gradient
{

template<class T>
LinearCombination<T, Vector> greenGauseGradientImpl
(
    MeshBase const& mesh,
    Index cellIdx,
    BoundaryConditionGetter<T> const& boundaries
)
{
    LinearCombination<T, Vector> gradient;

    for (Index faceIdx : mesh.getCellFaces(cellIdx))
    {
        Vector faceVector = mesh.getFaceVector(faceIdx);
        if (cellIdx != mesh.getFaceOwner(faceIdx))
        {
            faceVector *= -1;
        }
        // Order of multiplication is important
        gradient += valueOnFace(mesh, faceIdx, boundaries) * faceVector;
    }
    gradient /= mesh.getCellVolume(cellIdx);

    return gradient;
}


template<class T>
LinearCombination<T, Vector> leastSquareGradientImpl
(
    MeshBase const& mesh,
    Index cellIdx,
    BoundaryConditionGetter<T> const& boundaries
)
{
    Tensor systemMatrix = Tensor::Zero();
    LinearCombination<T, Vector> rhs;

    for (Index faceIdx : mesh.getCellFaces(cellIdx))
    {
        Vector radiusVector;
        LinearCombination<T, Scalar> fieldChange = {{-1, cellIdx}};
        if (mesh.isBoundaryFace(faceIdx))
        {
            radiusVector = mesh.getFaceCentroid(faceIdx) - mesh.getCellCentroid(cellIdx);
            fieldChange += valueOnFace(mesh, faceIdx, boundaries);
        }
        else
        {
            Index neighborIdx;
            {
                auto [idx1, idx2] = mesh.getFaceNeighbors(faceIdx);
                neighborIdx = idx1 + idx2 - cellIdx;
            }
            radiusVector = mesh.getCellCentroid(neighborIdx) - mesh.getCellCentroid(cellIdx);
            fieldChange += Term<Scalar>{1, neighborIdx};
        }

        Scalar weight = 1 / radiusVector.norm();
        systemMatrix += weight * outerProduct(radiusVector, radiusVector);
        rhs += Vector(weight * radiusVector) * fieldChange;
    }

    systemMatrix(2,2) = 1;
    Tensor systemMatrixInverse = systemMatrix.inverse();

    // Inplace solving system
    rhs.bias = innerProduct(rhs.bias, systemMatrixInverse);
    for (auto& [coeff, _] : rhs.terms)
    {
        coeff = innerProduct(coeff, systemMatrixInverse);
    }

    return rhs;
}


template<class T>
LinearCombination<T, Scalar> faceNormalGradient
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

    // general case
    auto [tmpIdx1, tmpIdx2] = mesh.getFaceNeighbors(faceIdx);
    Index cellToIdx = (tmpIdx1 == cellFromIdx ? tmpIdx2 : tmpIdx1);

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

    auto cellFromGradient = leastSquareGradientImpl(mesh, cellFromIdx, boundaries);
    auto cellToGradient = leastSquareGradientImpl(mesh, cellToIdx, boundaries);

    Vector unitDirection = Geometry::cellToCellUnitVector(mesh, cellFromIdx, cellToIdx);
    Scalar cellFromDistanceToFace = Geometry::distanceCellToFaceInDirection(mesh, cellFromIdx, faceIdx, unitDirection);

    Scalar cellFromFactor = cellFromDistanceToFace / distanceBetweenCells;
    auto faceGradient = cellFromFactor * cellFromGradient + (1 - cellFromFactor) * cellToGradient;

    Vector faceVector = mesh.getFaceVector(faceIdx).normalized();
    if (mesh.getFaceOwner(faceIdx) != cellFromIdx)
    {
        faceVector *= -1;
    }

    // Using over-relaxed approach
    Vector orthogonalComponent = faceVector.norm() / faceVector.dot(unitDirection) * unitDirection;
    Vector nonOrthogonalComponent = faceVector - orthogonalComponent;
    
    return 
    (
        finiteDifference * orthogonalComponent.norm() +
        faceGradient.dot(nonOrthogonalComponent)
    );
}


} // namespace Interpolation::Schemes::Gradient

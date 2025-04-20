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

    // For 2D case matrix is degenerate
    if (mesh.is2D())
    {
        systemMatrix(2,2) = 1;
    }
    Tensor systemMatrixInverse = systemMatrix.inverse();

    // Inplace solving system
    rhs.bias = innerProduct(rhs.bias, systemMatrixInverse);
    for (auto& [coeff, _] : rhs.terms)
    {
        coeff = innerProduct(coeff, systemMatrixInverse);
    }

    return rhs;
}

} // namespace Interpolation::Schemes::Gradient

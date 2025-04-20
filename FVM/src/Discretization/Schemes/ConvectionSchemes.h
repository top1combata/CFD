#pragma once

#include "GradientSchemes.h"
#include "GradientComputation.h"
#include "Discretization/LinearCombination.h"
#include "Mesh/Geometry.h"
#include <iostream>

namespace Interpolation::Schemes::Convection
{

template<class T>
LinearCombination<T, Scalar> centralDifferenceImpl
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    [[maybe_unused]] Field<Scalar> const& massFlow
)
{
    return valueOnFace(mesh, faceIdx, boundaries);
}


template<class T>
LinearCombination<T, Scalar> upwindImpl
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Field<Scalar> const& massFlow
)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        return valueOnFace(mesh, faceIdx, boundaries);
    }

    auto [ownerIdx, neighbourIdx] = mesh.getFaceNeighbors(faceIdx);
    if (massFlow(faceIdx) < 0)
    {
        std::swap(ownerIdx, neighbourIdx);
    }

    return LinearCombination<T>({{1, ownerIdx}});
}


template<class T>
LinearCombination<T, Scalar> downwindImpl
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Field<Scalar> const& massFlow
)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        return valueOnFace(mesh, faceIdx, boundaries);
    }

    auto [ownerIdx, neighbourIdx] = mesh.getFaceNeighbors(faceIdx);
    if (massFlow(faceIdx) < 0)
    {
        std::swap(ownerIdx, neighbourIdx);
    }

    return LinearCombination<T>({{1, neighbourIdx}});
}


template<class T>
LinearCombination<T, Scalar> frommImpl
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Field<Scalar> const& massFlow
)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        return valueOnFace(mesh, faceIdx, boundaries);
    }
    
    auto [ownerIdx, neighbourIdx] = mesh.getFaceNeighbors(faceIdx);
    if (massFlow(faceIdx) < 0)
    {
        std::swap(ownerIdx, neighbourIdx);
    }

    auto gradientScheme = 
    (
        mesh.useNonOrthogonalCorrection
        ? Gradient::LEAST_SQAURE
        : Gradient::GREEN_GAUSE
    );
    auto cellGradient = computeCellGradient(mesh, ownerIdx, boundaries, gradientScheme);
    Vector delta = mesh.getFaceCentroid(faceIdx) - mesh.getCellCentroid(ownerIdx);

    LinearCombination<T, Scalar> faceValue = {{1, ownerIdx}};
    faceValue += cellGradient.dot(delta);

    return faceValue + cellGradient.dot(delta);
}


template<class T>
LinearCombination<T, Scalar> souImpl
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    Field<Scalar> const& massFlow
)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        return valueOnFace(mesh, faceIdx, boundaries);
    }
    
    auto [ownerIdx, neighbourIdx] = mesh.getFaceNeighbors(faceIdx);
    if (massFlow(faceIdx) < 0)
    {
        std::swap(ownerIdx, neighbourIdx);
    }

    auto gradientScheme = 
    (
        mesh.useNonOrthogonalCorrection
        ? Gradient::LEAST_SQAURE
        : Gradient::GREEN_GAUSE
    );
    LinearCombination<T, Scalar> ownerValue = {{1, ownerIdx}};
    auto ownerCellGradient = computeCellGradient(mesh, ownerIdx, boundaries, gradientScheme);
    auto faceGradient = computeFaceGradient(mesh, faceIdx, boundaries, gradientScheme);
    Vector delta = mesh.getFaceCentroid(faceIdx) - mesh.getCellCentroid(ownerIdx);

    return ownerValue + Scalar(0.5) * (Scalar(2) * ownerCellGradient - faceGradient).dot(delta);
}


template<class T>
LinearCombination<T, Scalar> quickImpl
(
    MeshBase const& mesh,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries,
    [[maybe_unused]] Field<Scalar> const& massFlow
)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        return valueOnFace(mesh, faceIdx, boundaries);
    }
    
    auto [ownerIdx, neighbourIdx] = mesh.getFaceNeighbors(faceIdx);
    if (massFlow(faceIdx) < 0)
    {
        std::swap(ownerIdx, neighbourIdx);
    }

    auto gradientScheme = 
    (
        mesh.useNonOrthogonalCorrection
        ? Gradient::LEAST_SQAURE
        : Gradient::GREEN_GAUSE
    );
    LinearCombination<T, Scalar> ownerValue = {{1, ownerIdx}};
    auto ownerCellGradient = computeCellGradient(mesh, ownerIdx, boundaries, gradientScheme);
    auto faceGradient = computeFaceGradient(mesh, faceIdx, boundaries, gradientScheme);
    Vector delta = mesh.getFaceCentroid(faceIdx) - mesh.getCellCentroid(ownerIdx);

    return ownerValue + Scalar(0.5) * (ownerCellGradient + faceGradient).dot(delta);
}

} // nemespace Interpolation::Schemes::Convection

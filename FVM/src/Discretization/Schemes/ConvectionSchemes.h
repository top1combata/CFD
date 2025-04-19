#pragma once

#include "Mesh/MeshBase.h"
#include "Mesh/Geometry.h"
#include "Discretization/LinearCombination.h"
#include "BasicInterpolation.h"


namespace Interpolation::Schemes::Convection
{

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

} // nemespace Interpolation::Schemes::Convection

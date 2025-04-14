#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"
#include "Discretization/LinearCombination.h"

namespace Interpolation
{

template<class T>
LinearCombination<T, Vector> cellGradient
(
    MeshBase const& mesh, 
    Index cellIdx, 
    BoundaryConditionGetter<T> const& boundaries
);


inline Vector RhieChowVelocityOnFace
(
    MeshBase const& mesh,
    Index faceIdx,
    Field<Vector> const& U,
    Field<Scalar> const& p,
    Field<Vector> const& pGrad,
    Field<Scalar> const& VbyA,
    BoundaryConditionGetter<Vector> const& uBoundaries,
    BoundaryConditionGetter<Scalar> const& pBoundaries
);


template<class T>
LinearCombination<T, Scalar> valueOnFace
(
    MeshBase const& mesh, 
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries
);


template<class T>
LinearCombination<T, Scalar> faceNormalGradient
(
    MeshBase const& mesh, 
    Index cellFromIdx,
    Index faceIdx,
    BoundaryConditionGetter<T> const& boundaries
);


template<class T>
LinearCombination<T, Scalar> convectionFluxOverCell
(
    MeshBase const& mesh, 
    Index cellIdx, 
    BoundaryConditionGetter<T> const& boundaries, 
    Field<Scalar> const& massFlow
);


template<class T>
LinearCombination<T, Scalar> diffusionFluxOverCell
(
    MeshBase const& mesh, 
    Index cellIdx, 
    BoundaryConditionGetter<T> const& boundaries
);

} // namespace Interpolation

#include "Interpolation.hpp"

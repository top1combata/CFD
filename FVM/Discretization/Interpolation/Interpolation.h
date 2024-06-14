#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"
#include "Discretization/LinearCombination/LinearCombination.h"

namespace Interpolation
{

inline Vector cellGradient(MeshBase const& mesh, Index cellIdx, ScalarField const& field, BoundaryConditionGetter<Scalar> const& boundaries);

template<class T>
LinearCombination<T> valueOnFace(MeshBase const& mesh, Index faceIdx, BoundaryConditionGetter<T> const& boundaries);

template<class T>
LinearCombination<T> faceNormalGradient(MeshBase const& mesh, Index cellFromIdx, Index faceIdx, BoundaryConditionGetter<T> const& boundaries);

template<class T>
LinearCombination<T> convectionFluxOverCell(MeshBase const& mesh, Index cellIdx, BoundaryConditionGetter<T> const& boundaries, VectorField const& U);

template<class T>
LinearCombination<T> diffusionFluxOverCell(MeshBase const& mesh, Index cellIdx, BoundaryConditionGetter<T> const& boundaries);

} // namespace Interpolation

#include "Interpolation.hpp"
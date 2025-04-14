#pragma once

#include "Utils/Types.h"
#include "MeshBase.h"


namespace Geometry
{

Scalar distanceCellToFace(MeshBase const& mesh, Index cellIdx, Index faceIdx);

Scalar distanceCellToFaceInDirection(MeshBase const& mesh, Index cellIdx, Index faceIdx, Vector direction);

Scalar distanceCellToCell(MeshBase const& mesh, Index cellFromIdx, Index cellToIdx);

Vector cellToCellUnitVector(MeshBase const& mesh, Index cellFromIdx, Index cellToIdx);

} // namespace Geometry

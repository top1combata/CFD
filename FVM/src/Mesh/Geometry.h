#pragma once

#include "Utils/Types.h"
#include "MeshBase.h"


namespace Geometry
{

Scalar distanceCellToFace(MeshBase const& mesh, Index cellIdx, Index faceIdx);

Scalar distanceCellToCell(MeshBase const& mesh, Index cellFromIdx, Index cellToIdx);

}
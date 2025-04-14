#include "Geometry.h"


namespace Geometry
{

Scalar distanceCellToFace(MeshBase const& mesh, Index cellIdx, Index faceIdx)
{
    Vector normal = mesh.getFaceVector(faceIdx);
    return distanceCellToFaceInDirection(mesh, cellIdx, faceIdx, normal);
}

Scalar distanceCellToFaceInDirection(MeshBase const& mesh, Index cellIdx, Index faceIdx, Vector direction)
{
    direction.normalize();
    Vector cellCentroid = mesh.getCellCentroid(cellIdx);
    Vector faceCentroid = mesh.getFaceCentroid(faceIdx);
    Vector faceNormal = mesh.getFaceVector(faceIdx);

    Scalar dist = abs(faceNormal.dot(faceCentroid) - faceNormal.dot(cellCentroid)) / faceNormal.dot(direction);
    return dist;
}

Scalar distanceCellToCell(MeshBase const& mesh, Index cellFromIdx, Index cellToIdx)
{
    Vector radius = mesh.getCellCentroid(cellFromIdx) - mesh.getCellCentroid(cellToIdx);
    return radius.norm();
}

Vector cellToCellUnitVector(MeshBase const& mesh, Index cellFromIdx, Index cellToIdx)
{
    Vector radius = mesh.getCellCentroid(cellToIdx) - mesh.getCellCentroid(cellFromIdx);
    return radius.normalized();
}

} // namespace Geometry

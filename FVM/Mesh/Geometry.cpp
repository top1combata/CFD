#include "Geometry.h"


Scalar Geometry::distanceCellToFace(MeshBase const& mesh, Index cellIdx, Index faceIdx)
{
    Vector faceCentroid = mesh.getFaceCentroid(faceIdx);
    Vector cellCentroid = mesh.getCellCentroid(cellIdx);
    Vector normal = mesh.getFaceVector(faceIdx);

    Scalar dist = abs(cellCentroid.dot(normal) - faceCentroid.dot(normal)) / normal.norm();
    return dist;
}

Scalar Geometry::distanceCellToCell(MeshBase const& mesh, Index cellFromIdx, Index cellToIdx)
{
    Vector radius = mesh.getCellCentroid(cellFromIdx) - mesh.getCellCentroid(cellToIdx);
    return radius.norm();
}
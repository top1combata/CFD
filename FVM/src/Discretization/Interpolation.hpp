#include "Interpolation.h"

#include "Mesh/Geometry.h"


template<class T>
LinearCombination<T, Vector> Interpolation::cellGradient
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
            faceVector *= -1;
        gradient += valueOnFace(mesh, faceIdx, boundaries) * faceVector;
    }
    gradient /= mesh.getCellVolume(cellIdx);

    return gradient;
}


template<class T>
LinearCombination<T, Scalar> Interpolation::diffusionFluxOverCell
(
    MeshBase const& mesh, 
    Index cellIdx, 
    BoundaryConditionGetter<T> const& boundaries
)
{
    LinearCombination<T> flux;
    for (Index faceIdx : mesh.getCellFaces(cellIdx))
    {
        Vector faceVector = mesh.getFaceVector(faceIdx);
        if (mesh.getFaceOwner(faceIdx) != cellIdx)
            faceVector *= -1;
        flux += faceNormalGradient(mesh, cellIdx, faceIdx, boundaries) * faceVector.norm(); 
    }
    return flux;
}


template<class T>
LinearCombination<T, Scalar> Interpolation::convectionFluxOverCell
(
    MeshBase const& mesh, 
    Index cellIdx, 
    BoundaryConditionGetter<T> const& boundaries, 
    Field<Scalar> const& massFlow
)
{
    LinearCombination<T> flux;

    for (Index faceIdx : mesh.getCellFaces(cellIdx))
    {
        Scalar massFlux = massFlow(faceIdx);

        auto [ownerIdx, neighbourIdx] = mesh.getFaceNeighbors(faceIdx);
        
        LinearCombination<T> implicitVelocity = {{1,1}};
        implicitVelocity.terms[0].idx = (massFlux > 0 ? ownerIdx : neighbourIdx);
        if (mesh.isBoundaryFace(faceIdx))
            implicitVelocity = valueOnFace(mesh, faceIdx, boundaries);
        
        if (cellIdx != ownerIdx)
            massFlux *= -1;
        flux += massFlux * implicitVelocity;
    }

    return flux;
}


template<class T>
LinearCombination<T, Scalar> Interpolation::valueOnFace
(
    MeshBase const& mesh, 
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
        case FIXED_VALUE:
            return LinearCombination<T, Scalar>(boundaryValue);

        case FIXED_GRADIENT:
            Index cellIdx = mesh.getFaceNeighbors(faceIdx).front();
            Scalar dist = Geometry::distanceCellToFace(mesh, cellIdx, faceIdx);

            LinearCombination<T> result = {{1, cellIdx}};
            result += dist*boundaryValue;
            return  result;
        }

    }

    // general case
    auto [cellFromIdx, cellToIdx] = mesh.getFaceNeighbors(faceIdx);

    Scalar coeffFrom, coeffTo;

    if (mesh.useNonOrthogonalCorrection)
    {
        Scalar distanceBetweenCells = Geometry::distanceCellToCell(mesh, cellFromIdx, cellToIdx);
        Vector unitDirection = Geometry::cellToCellUnitVector(mesh, cellFromIdx, cellToIdx);
        Scalar distanceToFace = Geometry::distanceCellToFaceInDirection(mesh, cellFromIdx, faceIdx, unitDirection);
    
        coeffFrom = distanceToFace / distanceBetweenCells;
        coeffTo = 1 - coeffFrom;
    }
    else 
    {
        Scalar distFrom = Geometry::distanceCellToFace(mesh, cellFromIdx, faceIdx);
        Scalar distTo = Geometry::distanceCellToFace(mesh, cellToIdx, faceIdx);
    
        Scalar coeffFrom = 1 - distFrom/ (distFrom + distTo);
        Scalar coeffTo = 1 - coeffFrom;
    }
    return {{coeffFrom, cellFromIdx}, {coeffTo, cellToIdx}};
}


template<class T>
LinearCombination<T, Scalar> Interpolation::faceNormalGradient
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

    auto cellFromGradient = cellGradient(mesh, cellFromIdx, boundaries);
    auto cellToGradient = cellGradient(mesh, cellToIdx, boundaries);

    Vector unitDirection = Geometry::cellToCellUnitVector(mesh, cellFromIdx, cellToIdx);
    Scalar cellFromDistanceToFace = Geometry::distanceCellToFaceInDirection(mesh, cellFromIdx, faceIdx, unitDirection);

    Scalar cellFromFactor = cellFromDistanceToFace / distanceBetweenCells;
    auto faceGradient = cellFromFactor * cellFromGradient + (1 - cellFromFactor) * cellToGradient;

    Vector faceVector = mesh.getFaceVector(faceIdx).normalized();
    if (mesh.getFaceOwner(faceIdx) != cellFromIdx)
        faceVector *= -1;

    // Using over-relaxed approach
    Vector orthogonalComponent = faceVector.norm() / faceVector.dot(unitDirection) * unitDirection;
    Vector nonOrthogonalComponent = faceVector - orthogonalComponent;
    
    return 
    (
        finiteDifference * orthogonalComponent.norm() +
        faceGradient.dot(nonOrthogonalComponent)
    );
}


static Vector Interpolation::RhieChowVelocityOnFace
(
    MeshBase const& mesh,
    Index faceIdx,
    Field<Vector> const& U,
    Field<Scalar> const& p,
    Field<Vector> const& pGrad,
    Field<Scalar> const& VbyA,
    BoundaryConditionGetter<Vector> const& uBoundaries,
    BoundaryConditionGetter<Scalar> const& pBoundaries
)
{
    Vector faceVelocity = valueOnFace(mesh, faceIdx, uBoundaries).evaluate(U);
    
    if (mesh.isBoundaryFace(faceIdx))
        return faceVelocity;

    Scalar VbyA_f = valueOnFace(mesh, faceIdx, zeroGradGetter<Scalar>()).evaluate(VbyA);

    Index cellIdx = mesh.getFaceOwner(faceIdx);
    Scalar faceNormalGrad = faceNormalGradient(mesh, cellIdx, faceIdx, pBoundaries).evaluate(p);
    Vector avgFaceGradient = valueOnFace(mesh, faceIdx, zeroGradGetter<Vector>()).evaluate(pGrad);
    // correction
    Vector faceVector = mesh.getFaceVector(faceIdx);
    Vector unitNormal = faceVector.normalized();
    Vector velocityCorrection = -VbyA_f * (faceNormalGrad - avgFaceGradient.dot(unitNormal))*unitNormal;
    
    return faceVelocity + velocityCorrection;
}

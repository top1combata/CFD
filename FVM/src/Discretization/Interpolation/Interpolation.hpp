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
            
        gradient += faceVector * valueOnFace(mesh, faceIdx, boundaries);
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
    auto [cellIdx1, cellIdx2] = mesh.getFaceNeighbors(faceIdx);
    Scalar dist1 = Geometry::distanceCellToFace(mesh, cellIdx1, faceIdx);
    Scalar dist2 = Geometry::distanceCellToFace(mesh, cellIdx2, faceIdx);

    Scalar coeff1 = 1 - dist1 / (dist1 + dist2);
    Scalar coeff2 = 1 - coeff1;
    return {{coeff1, cellIdx1}, {coeff2, cellIdx2}};
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

    Scalar dist = Geometry::distanceCellToCell(mesh, cellFromIdx, cellToIdx);
    return {{1/dist, cellToIdx}, {-1/dist, cellFromIdx}};
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
    Vector unitNormal = faceVector / faceVector.norm();
    Vector velocityCorrection = -VbyA_f * (faceNormalGrad - avgFaceGradient.dot(unitNormal))*unitNormal;
    
    return faceVelocity + velocityCorrection;
}
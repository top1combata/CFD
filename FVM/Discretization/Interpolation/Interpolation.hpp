#include "Interpolation.h"
#include "Mesh/Geometry.h"


inline Vector Interpolation::cellGradient(MeshBase const& mesh, Index cellIdx, ScalarField const& field, BoundaryConditionGetter<Scalar> const& boundaries)
{
    Vector gradient(0,0,0);
    for (auto [faceVector, faceIdx] : mesh.getCellFaces(cellIdx))
    {
        gradient += faceVector * valueOnFace(mesh, faceIdx, boundaries).evaluate(field);
    }

    return 1/mesh.getCellVolume(cellIdx) * gradient;
}


template<class T>
LinearCombination<T> Interpolation::diffusionFluxOverCell(MeshBase const& mesh, Index cellIdx, BoundaryConditionGetter<T> const& boundaries)
{
    LinearCombination<T> flux;
    for (auto [faceVector, faceIdx] : mesh.getCellFaces(cellIdx))
    {
        flux += faceNormalGradient(mesh, cellIdx, faceIdx, boundaries) * faceVector.norm();;
    }
    return flux;
}


template<class T>
LinearCombination<T> Interpolation::convectionFluxOverCell(MeshBase const& mesh, Index cellIdx, BoundaryConditionGetter<T> const& boundaries, VectorField const& U)
{
    LinearCombination<T> flux;
    for (auto [faceVector, faceIdx] : mesh.getCellFaces(cellIdx))
    {
        auto velocityBoundaryGetter = [&mesh](Index idx)
        {
            return mesh.getFaceBoundary(idx).uBoundary;
        };
        Vector velocityOnFace = valueOnFace(mesh, faceIdx, velocityBoundaryGetter).evaluate(U);
        
        flux += velocityOnFace.dot(faceVector) * valueOnFace(mesh, faceIdx, boundaries);
    }

    return flux;
}


template<class T>
LinearCombination<T> Interpolation::valueOnFace(MeshBase const& mesh, Index faceIdx, BoundaryConditionGetter<T> const& boundaries)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        auto [boundaryValue, boundaryType] = boundaries(faceIdx);

        using enum BoundaryConditionType;

        switch (boundaryType)
        {
        case FIXED_VALUE:
            return boundaryValue;

        case FIXED_GRADIENT:
            Index cellIdx = mesh.getFaceNeighbours(faceIdx).front();
            Scalar dist = Geometry::distanceCellToFace(mesh, cellIdx, faceIdx);

            LinearCombination<T> result;
            result += dist*boundaryValue;
            result += {{1, cellIdx}};
            return  result;
        }

    }

    // general case
    auto [cellIdx1, cellIdx2] = mesh.getFaceNeighbours(faceIdx);
    Scalar dist1 = Geometry::distanceCellToFace(mesh, cellIdx1, faceIdx);
    Scalar dist2 = Geometry::distanceCellToFace(mesh, cellIdx2, faceIdx);

    Scalar coeff1 = 1 - dist1 / (dist1 + dist2);
    Scalar coeff2 = 1 - coeff1;
    return {{coeff1, cellIdx1}, {coeff2, cellIdx2}};
}


template<class T>
LinearCombination<T> Interpolation::faceNormalGradient(MeshBase const& mesh, Index cellFromIdx, Index faceIdx, BoundaryConditionGetter<T> const& boundaries)
{
    if (mesh.isBoundaryFace(faceIdx))
    {
        auto [boundaryValue, boundaryType] = boundaries(faceIdx);

        using enum BoundaryConditionType;

        switch (boundaryType)
        {
        case FIXED_GRADIENT:
            return boundaryValue;
        
        case FIXED_VALUE:
            Scalar dist = Geometry::distanceCellToFace(mesh, cellFromIdx, faceIdx);

            LinearCombination<T> result;
            result += boundaryValue / dist;
            result -= {{1/dist, cellFromIdx}};
            return result;
        }
    }

    // general case
    auto [tmpIdx1, tmpIdx2] = mesh.getFaceNeighbours(faceIdx);
    Index cellToIdx = (tmpIdx1 == cellFromIdx ? tmpIdx2 : tmpIdx1);

    Scalar dist = Geometry::distanceCellToCell(mesh, cellFromIdx, cellToIdx);
    return {{1/dist, cellToIdx}, {-1/dist, cellFromIdx}};
}
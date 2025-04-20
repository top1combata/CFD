#pragma once

#include "Schemes/InterpolationSchemes.h"
#include "Schemes/BasicInterpolation.h"
#include "Schemes/GradientSchemes.h"
#include "Schemes/GradientComputation.h"
#include "Schemes/ConvectionSchemes.h"


namespace Interpolation
{

template<class T>
LinearCombination<T, Scalar> computeConvectionFluxOverCell
(
    MeshBase const& mesh, 
    Index cellIdx, 
    BoundaryConditionGetter<T> const& boundaries, 
    Field<Scalar> const& massFlow,
    Schemes::Convection::Type schemeType
)
{
    LinearCombination<T> convectiveFlux;

    decltype(&Schemes::Convection::upwindImpl<T>) schemeImpl;
    switch (schemeType)
    {
        case Schemes::Convection::CENTRAL_DIFFERENCE:
            schemeImpl = &Schemes::Convection::centralDifferenceImpl;
            break;

        case Schemes::Convection::UPWIND:
            schemeImpl = &Schemes::Convection::upwindImpl;
            break;

        case Schemes::Convection::DOWNWIND:
            schemeImpl = &Schemes::Convection::downwindImpl;
            break;

        case Schemes::Convection::FROMM:
            schemeImpl = &Schemes::Convection::frommImpl;
            break;

        case Schemes::Convection::SOU:
            schemeImpl = &Schemes::Convection::souImpl;
            break;

        case Schemes::Convection::QUICK:
            schemeImpl = &Schemes::Convection::quickImpl;
            break;
    }

    for (Index faceIdx : mesh.getCellFaces(cellIdx))
    {
        Scalar massFlux = massFlow(faceIdx);

        if (cellIdx != mesh.getFaceOwner(faceIdx))
        {
            massFlux *= -1;
        }

        convectiveFlux += 
        (
            massFlux *
            (*schemeImpl)
            (
                mesh, faceIdx, boundaries, massFlow
            )
        );
    }

    return convectiveFlux;
}


template<class T>
LinearCombination<T, Scalar> computeDiffusionFluxOverCell
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
        {
            faceVector *= -1;
        }
        flux += faceVector.norm() *
        computeFaceNormalGradient
        (
            mesh, cellIdx, faceIdx, boundaries
        ); 
    }
    return flux;
}


inline Vector computeRhieChowVelocityOnFace
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
    {
        return faceVelocity;
    }

    Scalar VbyA_f = valueOnFace(mesh, faceIdx, zeroGradGetter<Scalar>()).evaluate(VbyA);

    Index cellIdx = mesh.getFaceOwner(faceIdx);
    Scalar faceNormalGrad = computeFaceNormalGradient(mesh, cellIdx, faceIdx, pBoundaries).evaluate(p);
    Vector avgFaceGradient = valueOnFace(mesh, faceIdx, zeroGradGetter<Vector>()).evaluate(pGrad);
    // correction
    Vector faceVector = mesh.getFaceVector(faceIdx);
    Vector unitNormal = faceVector.normalized();
    Vector velocityCorrection = -VbyA_f * (faceNormalGrad - avgFaceGradient.dot(unitNormal))*unitNormal;
    
    return faceVelocity + velocityCorrection;
}

} // namespace Interpolation

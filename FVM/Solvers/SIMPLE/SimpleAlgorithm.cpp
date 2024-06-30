#include "SimpleAlgorithm.h"
#include "Discretization/Interpolation/Interpolation.h"
#include "Config/PhysicalProperties.h"
#include "Config/SolverControl.h"
#include "Eigen/Dense"


SimpleAlgorithm::SimpleAlgorithm(MeshBase const& mesh) : m_mesh(mesh){}

#include <iostream>

void SimpleAlgorithm::solve()
{
    initFields();

    for (Index iterNum = 1; iterNum <= Config::maxIterations && !converged(); iterNum++)
    {
        computePressureGradient();
        solveMomentum();
        updateMassFluxes();
        correctPressure();

        std::cout << "Iteration #" << iterNum << 
            "\n\tU residual : " << m_U_residual << 
            "\n\tp residual : " << m_p_residual << 
            "\n\n";
    }
}



void SimpleAlgorithm::initFields()
{
    Index totalCells = m_mesh.getCellsAmount();
    Index totalFaces = m_mesh.getFacesAmount();

    m_p           = ScalarField::Zero(totalCells, 1);
    m_U           = VectorField::Zero(totalCells, 3);
    m_mass_fluxes = ScalarField::Zero(totalFaces, 1);
    m_p_grad      = VectorField::Zero(totalCells, 3);
    m_Au          = Matrix::Zero(totalCells, totalCells);
    m_VbyA        = ScalarField::Zero(totalCells, 1);

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Vector faceVelocity = Interpolation::valueOnFace(m_mesh, faceIdx, uBoundaries()).evaluate(m_U);
            m_mass_fluxes(faceIdx) = faceVelocity.dot(faceVector) * Config::density;
        }
    }
}


void SimpleAlgorithm::computePressureGradient()
{
    Index totalCells = m_mesh.getCellsAmount();

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_p_grad.row(cellIdx) = Interpolation::cellGradient(m_mesh, cellIdx, m_p, pBoundaries()).transpose();
    }
}


void SimpleAlgorithm::solveMomentum()
{
    Index totalCells = m_mesh.getCellsAmount();

    m_Au.setZero();
    VectorField uSource = VectorField::Zero(totalCells, 3);

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        LinearCombination<Vector> convection = Interpolation::convectionFluxOverCell(m_mesh, cellIdx, uBoundaries(), m_mass_fluxes);
        LinearCombination<Vector> diffusion  = Interpolation::diffusionFluxOverCell(m_mesh, cellIdx, uBoundaries()) * Config::viscosity;
        Vector pressureGradient = getFieldValue(m_p_grad, cellIdx) * m_mesh.getCellVolume(cellIdx);

        LinearCombination<Vector> Ueqn;
        Ueqn += convection;
        Ueqn -= diffusion;
        Ueqn += pressureGradient;   

        for (auto [coeff, idx] : Ueqn.terms)
            m_Au(cellIdx, idx) += coeff;

        uSource.row(cellIdx) =  -Ueqn.bias.transpose();
    }

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_VbyA(cellIdx) = m_mesh.getCellVolume(cellIdx) / m_Au(cellIdx, cellIdx);
    }

    // implicit relaxation
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        uSource.row(cellIdx) += (1/Config::uRelax - 1) * m_Au(cellIdx, cellIdx) * m_U.row(cellIdx);
        m_Au(cellIdx, cellIdx) /= Config::uRelax;
    }

    m_U = m_Au.partialPivLu().solve(uSource);
}


void SimpleAlgorithm::updateMassFluxes()
{
    Index totalCells = m_mesh.getCellsAmount();
    
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
            
            Vector faceVelocity = Interpolation::valueOnFace(m_mesh, faceIdx, uBoundaries()).evaluate(m_U);

            if (m_mesh.isBoundaryFace(faceIdx) == false)
            {
                Scalar faceNormalGradient = Interpolation::faceNormalGradient(m_mesh, cellIdx, faceIdx, pBoundaries()).evaluate(m_p);
                Vector avgFaceGradient = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Vector>()).evaluate(m_p_grad);
                // correction
                Vector unitNormal = faceVector / faceVector.norm();
                Vector velocityCorrection = -VbyAf * (faceNormalGradient - avgFaceGradient.dot(unitNormal))*unitNormal;
                
                faceVelocity += velocityCorrection;
            }
            
            m_mass_fluxes(faceIdx) = faceVelocity.dot(faceVector) * Config::density;
        }
    }
}


void SimpleAlgorithm::correctPressure()
{
    Index totalCells = m_mesh.getCellsAmount();
    Index totalFaces = m_mesh.getFacesAmount();


    BoundaryConditionGetter<Scalar> pCorrBoundaries = [pBound = pBoundaries()] (Index faceIdx)
    {
        auto boundaries = pBound(faceIdx);
        boundaries.value = 0;
        return boundaries;
    };

    Matrix pA = Matrix::Zero(totalCells, totalCells);
    ScalarField pSource = ScalarField::Zero(totalCells);

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {               
        LinearCombination<Scalar> diffusiveFlux;
        Scalar massFlow = 0;
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
            
            diffusiveFlux += VbyAf *Config::density * faceVector.norm() * Interpolation::faceNormalGradient(m_mesh, cellIdx, faceIdx, pCorrBoundaries);

            massFlow += m_mass_fluxes(faceIdx);
        }

        LinearCombination<Scalar> Peqn;
        Peqn -= massFlow;
        Peqn += diffusiveFlux;

        for (auto [coeff, idx] : Peqn.terms)
            pA(cellIdx, idx) += coeff;
        
        pSource(cellIdx) = -Peqn.bias;
    }

    ScalarField pCorrection = pA.partialPivLu().solve(pSource);


    VectorField uCorrection(totalCells, 3);
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
        uCorrection.row(cellIdx) = -m_VbyA(cellIdx) * Interpolation::cellGradient(m_mesh, cellIdx, pCorrection, pCorrBoundaries).transpose();

    ScalarField massFluxesCorrection(totalFaces, 1);
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
            massFluxesCorrection(faceIdx) = -Config::density * VbyAf * Interpolation::faceNormalGradient(m_mesh, cellIdx, faceIdx, pCorrBoundaries).evaluate(pCorrection);
        }
    }

    m_U_residual = relativeResidual(m_U, uCorrection);
    m_p_residual = relativeResidual(m_p, pCorrection);

    m_U           += uCorrection * Config::pRelax;
    m_p           += pCorrection * Config::pRelax;
    m_mass_fluxes += massFluxesCorrection * Config::pRelax;
}


bool SimpleAlgorithm::converged()
{
    return  m_U_residual < Config::uTolerance && m_p_residual < Config::pTolerance
        || std::isnan(m_U_residual) && std::isnan(m_p_residual);
}


Scalar SimpleAlgorithm::relativeResidual(Matrix const& field, Matrix const& correction)
{
    Scalar fieldMax = std::max(field.maxCoeff(), -field.minCoeff());
    Scalar corrMax  = std::max(correction.maxCoeff(), -correction.minCoeff());
    
    return corrMax / fieldMax;
}


BoundaryConditionGetter<Vector> SimpleAlgorithm::uBoundaries()
{
    return [&mesh = m_mesh](Index faceIdx)
    {
        return mesh.getFaceBoundary(faceIdx).uBoundary;
    }; 
}


BoundaryConditionGetter<Scalar> SimpleAlgorithm::pBoundaries()
{
    return [&mesh = m_mesh](Index faceIdx)
    {
        return mesh.getFaceBoundary(faceIdx).pBoundary;
    }; 
}
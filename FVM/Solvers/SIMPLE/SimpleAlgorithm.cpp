#include "SimpleAlgorithm.h"
#include "Discretization/Interpolation/Interpolation.h"
#include "Config/PhysicalProperties.h"
#include "Config/SolverControl.h"
#include "Utils/MatrixSolver.h"
#include <iostream>
#include <format>


SimpleAlgorithm::SimpleAlgorithm(MeshBase const& mesh) : m_mesh(mesh){}


void SimpleAlgorithm::solve()
{
    initFields();

    // SIMPLE loop
    m_timers["total"].start();
    for (Index iterNum = 1; iterNum <= Config::maxIterations && !converged(); iterNum++)
    {
        computePressureGradient();
        solveMomentum();
        computeMassFluxes();
        correctPressure();

        std::cout << "Iteration #" << iterNum << 
            "\n\tU residual : " << m_U_residual << 
            "\n\tp residual : " << m_p_residual << 
            "\n\n";
    }
    m_timers["total"].stop();

    // how much time elapsed
    for (auto&& [eventName, eventTimer] : m_timers)
        std::cout << std::format("{} seconds for {}\n\n", eventTimer.getElapsedTime(), eventName);
}


bool SimpleAlgorithm::converged()
{
    return  (m_U_residual < Config::uTolerance && m_p_residual < Config::pTolerance)
        || (std::isnan(m_U_residual) && std::isnan(m_p_residual));
}


void SimpleAlgorithm::solveMomentum()
{
    generateMomentumSystem();

    relaxSystem(m_U_matrix, m_U_source, m_U, Config::uRelax);

    // the field V/A should be treated after under relaxation
    computeVbyA();

    // solving for new velocity field
    m_timers["solving linear systems"].start();
    m_U = solveSystem(m_U_matrix, m_U_source);
    m_timers["solving linear systems"].stop();
}


void SimpleAlgorithm::correctPressure()
{
    generatePressureCorrectionSystem();

    m_timers["solving linear systems"].start();
    ScalarField pCorrection = solveSystem(m_p_matrix, m_p_source);
    m_timers["solving linear systems"].stop();

    // idk how, but explicit under relaxtion gives faster convergence than implicit
    pCorrection *= Config::pRelax;

    VectorField uCorrection = getVelocityCorrection(pCorrection);
    ScalarField massFluxesCorrection = getMassFluxesCorrection(pCorrection);

    m_U_residual = relativeResidual(m_U, uCorrection);
    m_p_residual = relativeResidual(m_p, pCorrection);

    m_U           += uCorrection;
    m_p           += pCorrection;
    m_mass_fluxes += massFluxesCorrection;
}


VectorField& SimpleAlgorithm::getU()
{
    return m_U;
}


ScalarField& SimpleAlgorithm::getP()
{
    return m_p;
}


void SimpleAlgorithm::initFields()
{
    Index totalCells = m_mesh.getCellAmount();
    Index totalFaces = m_mesh.getFaceAmount();

    m_p           = ScalarField::Zero(totalCells, 1);
    m_U           = VectorField::Zero(totalCells, 3);
    m_mass_fluxes = ScalarField::Zero(totalFaces, 1);
    m_p_grad      = VectorField::Zero(totalCells, 3);
    m_U_matrix    = Matrix::Zero(totalCells, totalCells);
    m_U_source    = Matrix::Zero(totalCells, 3);
    m_p_matrix    = Matrix::Zero(totalCells, totalCells);
    m_p_source    = Matrix::Zero(totalCells, 1);
    m_VbyA        = ScalarField::Zero(totalCells, 1);

    m_timers.clear();

    // init mass fluxes
    // can't be just zero because of boundary values
    auto uBoundaries = getVelocityBoundaries();
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Vector faceVelocity = Interpolation::valueOnFace(m_mesh, faceIdx, uBoundaries).evaluate(m_U);
            m_mass_fluxes(faceIdx) = faceVelocity.dot(faceVector) * Config::density;
        }
    }
}


void SimpleAlgorithm::computePressureGradient()
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pBoundaries = getPressureBoundaries();

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_p_grad.row(cellIdx) = Interpolation::cellGradient(m_mesh, cellIdx, m_p, pBoundaries).transpose();
    }
    m_timers["explicit field computation"].stop();
}


void SimpleAlgorithm::computeMassFluxes()
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pBoundaries = getPressureBoundaries();
    auto uBoundaries = getVelocityBoundaries();
    
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
            
            Vector faceVelocity = Interpolation::valueOnFace(m_mesh, faceIdx, uBoundaries).evaluate(m_U);

            if (m_mesh.isBoundaryFace(faceIdx) == false)
            {
                Scalar faceNormalGradient = Interpolation::faceNormalGradient(m_mesh, cellIdx, faceIdx, pBoundaries).evaluate(m_p);
                Vector avgFaceGradient = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Vector>()).evaluate(m_p_grad);
                // correction
                Vector unitNormal = faceVector / faceVector.norm();
                Vector velocityCorrection = -VbyAf * (faceNormalGradient - avgFaceGradient.dot(unitNormal))*unitNormal;
                
                faceVelocity += velocityCorrection;
            }
            
            m_mass_fluxes(faceIdx) = faceVelocity.dot(faceVector) * Config::density;
        }
    }
    m_timers["explicit field computation"].stop();
}


void SimpleAlgorithm::generateMomentumSystem()
{
    m_timers["generating linear systems"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto uBoundaries = getVelocityBoundaries();

    m_U_matrix.setZero();
    m_U_source.setZero();

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        LinearCombination<Vector> convection = Interpolation::convectionFluxOverCell(m_mesh, cellIdx, uBoundaries, m_mass_fluxes);
        LinearCombination<Vector> diffusion  = Interpolation::diffusionFluxOverCell(m_mesh, cellIdx, uBoundaries) * Config::viscosity;
        Vector pressureGradient = getFieldValue(m_p_grad, cellIdx) * m_mesh.getCellVolume(cellIdx);

        LinearCombination<Vector> Ueqn;
        Ueqn += convection;
        Ueqn -= diffusion;
        Ueqn += pressureGradient;   

        for (auto [coeff, idx] : Ueqn.terms)
            m_U_matrix(cellIdx, idx) += coeff;

        m_U_source.row(cellIdx) =  -Ueqn.bias.transpose();
    }
    m_timers["generating linear systems"].stop();
}


void SimpleAlgorithm::generatePressureCorrectionSystem()
{
    m_timers["generating linear systems"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    m_p_matrix.setZero();
    m_p_source.setZero();

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
            m_p_matrix(cellIdx, idx) += coeff;
        
        m_p_source(cellIdx) = -Peqn.bias;
    }
    m_timers["generating linear systems"].stop();
}


void SimpleAlgorithm::computeVbyA()
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_VbyA(cellIdx) = m_mesh.getCellVolume(cellIdx) / m_U_matrix(cellIdx, cellIdx);
    }
    m_timers["explicit field computation"].stop();
}


VectorField SimpleAlgorithm::getVelocityCorrection(ScalarField const& pCorrection)
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    VectorField uCorrection(totalCells, 3);

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
        uCorrection.row(cellIdx) = -m_VbyA(cellIdx) * Interpolation::cellGradient(m_mesh, cellIdx, pCorrection, pCorrBoundaries).transpose();
    
    m_timers["explicit field computation"].stop();
    return uCorrection;
}


ScalarField SimpleAlgorithm::getMassFluxesCorrection(ScalarField const& pCorrection)
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    Index totalFaces = m_mesh.getFaceAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    ScalarField massFluxesCorrection(totalFaces, 1);

    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (auto [faceVector, faceIdx] : m_mesh.getCellFaces(cellIdx))
        {
            Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
            massFluxesCorrection(faceIdx) = -Config::density * VbyAf * Interpolation::faceNormalGradient(m_mesh, cellIdx, faceIdx, pCorrBoundaries).evaluate(pCorrection);
        }
    }
    m_timers["explicit field computation"].stop();
    return massFluxesCorrection; 
}


Scalar SimpleAlgorithm::relativeResidual(Matrix const& field, Matrix const& correction)
{
    m_timers["computing residuals"].start();
    Scalar fieldMax = std::max(field.maxCoeff(), -field.minCoeff());
    Scalar corrMax  = std::max(correction.maxCoeff(), -correction.minCoeff());
    m_timers["computing residuals"].stop();
    
    return corrMax / fieldMax;
}


BoundaryConditionGetter<Vector> SimpleAlgorithm::getVelocityBoundaries()
{
    return [&mesh = m_mesh](Index faceIdx)
    {
        return mesh.getFaceBoundary(faceIdx).uBoundary;
    }; 
}


BoundaryConditionGetter<Scalar> SimpleAlgorithm::getPressureBoundaries()
{
    return [&mesh = m_mesh](Index faceIdx)
    {
        return mesh.getFaceBoundary(faceIdx).pBoundary;
    }; 
}


BoundaryConditionGetter<Scalar> SimpleAlgorithm::getPressureCorrectionBoundaries()
{
    return [pBound = getPressureBoundaries()] (Index faceIdx)
    {
        auto boundaries = pBound(faceIdx);
        boundaries.value = 0;
        return boundaries;
    };
}
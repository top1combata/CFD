#include "SimpleAlgorithm.h"
#include "Discretization/Interpolation/Interpolation.h"
#include "Config/PhysicalProperties.h"
#include "Config/SolverControl.h"
#include "Utils/MatrixSolver.h"
#include <iostream>
#include <format>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif


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
    m_timers["generating linear systems"].start();
    generateMomentumSystem();
    m_timers["generating linear systems"].stop();

    relaxSystem(m_U_matrix, m_U_source, m_U, Config::uRelax);

    // the field V/A should be treated after under relaxation
    computeVbyA();

    // solving for new velocity field
    m_timers["solving linear systems"].start();
    m_U = solveSystem(m_U_matrix, m_U_source, m_U);
    m_timers["solving linear systems"].stop();
}


void SimpleAlgorithm::correctPressure()
{
    m_timers["generating linear systems"].start();
    generatePressureCorrectionSystem();
    m_timers["generating linear systems"].stop();

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

    // initial guesses
    m_p           = ScalarField::Zero(totalCells, 1);
    m_U           = VectorField::Zero(totalCells, 3);

    m_mass_fluxes = ScalarField(totalFaces, 1);
    m_p_grad      = VectorField(totalCells, 3);
    m_U_matrix    = SparseMatrix(totalCells, totalCells);
    m_U_source    = Matrix(totalCells, 3);
    m_p_matrix    = SparseMatrix(totalCells, totalCells);
    m_p_source    = Matrix(totalCells, 1);
    m_VbyA        = ScalarField(totalCells, 1);

    m_timers.clear();

    // init mass fluxes
    // can't be just zero because of boundary values
    auto uBoundaries = getVelocityBoundaries();
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (Index faceIdx : m_mesh.getCellFaces(cellIdx))
        {
            if (m_mesh.isBoundaryFace(faceIdx) && m_mesh.getFaceBoundary(faceIdx).uBoundary.type == BoundaryConditionType::FIXED_VALUE)
                m_mass_fluxes(faceIdx) = m_mesh.getFaceBoundary(faceIdx).uBoundary.value.dot(m_mesh.getFaceVector(faceIdx)) * Config::density;
            else
                m_mass_fluxes(faceIdx) = 0;
        }
    }
}


void SimpleAlgorithm::computePressureGradient()
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pBoundaries = getPressureBoundaries();

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_p_grad.row(cellIdx) = Interpolation::cellGradient(m_mesh, cellIdx, m_p, pBoundaries).transpose();
    }
    m_timers["explicit field computation"].stop();
}


void SimpleAlgorithm::computeMassFluxes()
{
    m_timers["explicit field computation"].start();
    Index totalFaces = m_mesh.getFaceAmount();
    auto pBoundaries = getPressureBoundaries();
    auto uBoundaries = getVelocityBoundaries();
    
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Vector faceVector = m_mesh.getFaceVector(faceIdx);
        Vector faceVelocity = Interpolation::RhieChowVelocityOnFace(m_mesh, faceIdx, m_U, m_p, m_p_grad, m_VbyA, uBoundaries, pBoundaries);
        m_mass_fluxes(faceIdx) = faceVelocity.dot(faceVector) * Config::density;
    }

    m_timers["explicit field computation"].stop();
}


template<class T, class Rhs>
void generateSparseSystemImpl(SparseMatrix&, Rhs&, Index size, std::function<LinearCombination<T>(Index)> const&);


void SimpleAlgorithm::generateMomentumSystem()
{
    std::function<LinearCombination<Vector>(Index)> uEqnGetter = 
    [this, uBoundaries = getVelocityBoundaries()](Index cellIdx)
    {
        LinearCombination<Vector> convection = Interpolation::convectionFluxOverCell(m_mesh, cellIdx, uBoundaries, m_mass_fluxes);
        LinearCombination<Vector> diffusion  = Interpolation::diffusionFluxOverCell(m_mesh, cellIdx, uBoundaries) * Config::viscosity;
        Vector pressureGradient = getFieldValue(m_p_grad, cellIdx) * m_mesh.getCellVolume(cellIdx);

        LinearCombination<Vector> uEqn;
        uEqn += convection;
        uEqn -= diffusion;
        uEqn += pressureGradient;

        return uEqn;
    };

    generateSparseSystemImpl(m_U_matrix, m_U_source, m_mesh.getCellAmount(), uEqnGetter);
}


void SimpleAlgorithm::generatePressureCorrectionSystem()
{
    std::function<LinearCombination<Scalar>(Index)> pCorrEqnGetter = 
    [this, pCorrBoundaries = getPressureCorrectionBoundaries()](Index cellIdx)
    {
        LinearCombination<Scalar> diffusiveFlux;
        Scalar massFlow = 0;
        for (Index faceIdx : m_mesh.getCellFaces(cellIdx))
        {
            Vector faceVector = m_mesh.getFaceVector(faceIdx);
            Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
            
            diffusiveFlux += VbyAf *Config::density * faceVector.norm() * Interpolation::faceNormalGradient(m_mesh, cellIdx, faceIdx, pCorrBoundaries);

            Scalar massFlux = m_mass_fluxes(faceIdx);
            if (cellIdx != m_mesh.getFaceOwner(faceIdx))
                massFlux *= -1;

            massFlow += massFlux;
        }

        LinearCombination<Scalar> pCorrEqn;
        pCorrEqn -= massFlow;
        pCorrEqn += diffusiveFlux;

        return pCorrEqn;
    };
    
    generateSparseSystemImpl(m_p_matrix, m_p_source, m_mesh.getCellAmount(), pCorrEqnGetter);
}


void SimpleAlgorithm::computeVbyA()
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto diagonal = m_U_matrix.diagonal();

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_VbyA(cellIdx) = m_mesh.getCellVolume(cellIdx) / diagonal(cellIdx);
    }
    m_timers["explicit field computation"].stop();
}


VectorField SimpleAlgorithm::getVelocityCorrection(ScalarField const& pCorrection)
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    VectorField uCorrection(totalCells, 3);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
        uCorrection.row(cellIdx) = -m_VbyA(cellIdx) * Interpolation::cellGradient(m_mesh, cellIdx, pCorrection, pCorrBoundaries).transpose();
    
    m_timers["explicit field computation"].stop();
    return uCorrection;
}


ScalarField SimpleAlgorithm::getMassFluxesCorrection(ScalarField const& pCorrection)
{
    m_timers["explicit field computation"].start();
    Index totalFaces = m_mesh.getFaceAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    ScalarField massFluxesCorrection(totalFaces, 1);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Index ownerIdx = m_mesh.getFaceOwner(faceIdx);
        Scalar VbyAf = Interpolation::valueOnFace(m_mesh, faceIdx, zeroGrad<Scalar>()).evaluate(m_VbyA);
        massFluxesCorrection(faceIdx) = -Config::density * VbyAf * Interpolation::faceNormalGradient(m_mesh, ownerIdx, faceIdx, pCorrBoundaries).evaluate(pCorrection);
    }
    m_timers["explicit field computation"].stop();
    return massFluxesCorrection; 
}


template<class T>
auto transpose(T vec)
{
    return vec.transpose();
}

template<>
auto transpose<Scalar>(Scalar value)
{
    return Eigen::Matrix<Scalar, 1, 1>(value);
}

template<class T, class Rhs>
void generateSparseSystemImpl(SparseMatrix& A, Rhs& rhs, Index size, std::function<LinearCombination<T>(Index)> const& eqnGetter)
{
    A.setZero();
    using Triplet = Eigen::Triplet<Scalar>;
    List<Triplet> triplets;

    auto cmp = [](Term lhs, Term rhs)
    {
        return lhs.idx < rhs.idx;
    };

#ifdef _OPENMP
    List<Triplet> threadTriplets;
    Index numThreads;
    #pragma omp parallel 
    #pragma omp single 
    numThreads = omp_get_num_threads();

    List<Index> sizes(numThreads);
    List<Index> prefix(numThreads+1, 0);

    #pragma omp parallel private(threadTriplets)
    {
        Index threadIdx = omp_get_thread_num();
        Index workRangeSize  = size / numThreads + bool(threadIdx < size % numThreads);
        Index workRangeStart = size / numThreads * threadIdx + std::min(threadIdx, size % numThreads);

        for (Index eqnIdx = workRangeStart; eqnIdx < workRangeStart + workRangeSize; eqnIdx++)
        {
            auto eqn = eqnGetter(eqnIdx);

            std::sort(eqn.terms.begin(), eqn.terms.end(), cmp);
            for (auto [coeff, varIdx] : eqn.terms)
                threadTriplets.emplace_back(eqnIdx, varIdx, coeff);

            rhs.row(eqnIdx) = -transpose(eqn.bias);
        }
        #pragma omp barrier
        sizes[threadIdx] = threadTriplets.size();
        #pragma omp barrier
        
        #pragma omp single
        {
            for (Index idx = 1; idx <= numThreads; idx++)
                prefix[idx] = prefix[idx-1] + sizes[idx-1];
            
            triplets.resize(prefix.back());
        }

        for (Index idx = 0; idx < sizes[threadIdx]; idx++)
            triplets[idx + prefix[threadIdx]] = threadTriplets[idx];
    }
#else
    for (Index eqnIdx = 0; eqnIdx < size; eqnIdx++)
    {
        auto eqn = eqnGetter(eqnIdx);

        std::sort(eqn.terms.begin(), eqn.terms.end(), cmp);
        for (auto [coeff, varIdx] : eqn.terms)
            triplets.emplace_back(eqnIdx, varIdx, coeff);

        rhs.row(eqnIdx) = -transpose(eqn.bias);
    }
#endif
    A.setFromSortedTriplets(triplets.begin(), triplets.end());
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
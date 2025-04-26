#include "SimpleAlgorithm.h"
#include "Discretization/Interpolation.h"
#include "Config/Config.h"
#include "Utils/MatrixSolver.h"
#include <iostream>
#include <format>
#include <algorithm>
#include <cassert>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif


SimpleAlgorithm::SimpleAlgorithm(MeshBase const& mesh) 
    : SolverBase(mesh)
    , gradientScheme(Config::gradientScheme)
    , convectionScheme(Config::convectionScheme)
{}


void SimpleAlgorithm::solve()
{
    initFields();

    // SIMPLE loop
    m_timers["total"].start();
    for 
    (
        Index iterNum = 1; 
        iterNum <= Config::maxIterations && !(converged() || diverged()); 
        iterNum++
    )
    {
        computePressureGradient();
        solveMomentum();
        computeMassFluxes();
        correctPressure();

        std::cout << "Iteration #" << iterNum << 
            "\n\tp residual : " << m_pressureResidual << 
            "\n\n";
    }

    m_velocity.push_back(m_currentVelocity);
    m_pressure.push_back(m_currentPressure);
    m_converged = converged();

    m_timers["total"].stop();

    // How much time elapsed
    for (auto&& [eventName, eventTimer] : m_timers)
    {
        std::cout << std::format("{} seconds for {}\n\n", eventTimer.getElapsedTime(), eventName);
    }
}


bool SimpleAlgorithm::converged() const
{
    return  
    (
        m_pressureResidual && m_pressureResidual < Config::pTolerance
    );
}


bool SimpleAlgorithm::diverged() const
{
    return  
    (
        std::isnan(m_pressureResidual) || m_pressureResidual == 0
    );
}


void SimpleAlgorithm::solveMomentum()
{
    m_timers["generating linear systems"].start();
    generateMomentumSystem();
    m_timers["generating linear systems"].stop();

    relaxSystem(m_momentumSystemMatrix, m_momentumSystemSource, m_currentVelocity, Config::uRelax);

    // Field V/A should be treated after under relaxation
    computeVbyA();

    // Solving for new velocity field
    m_timers["solving linear systems"].start();
    auto sol = solveSystem(m_momentumSystemMatrix, m_momentumSystemSource);
    for (Index cellIdx = 0; cellIdx < m_mesh.getCellAmount(); cellIdx++)
    {
        m_currentVelocity(cellIdx) = sol.row(cellIdx).transpose();
    }
    m_timers["solving linear systems"].stop();
}


void SimpleAlgorithm::correctPressure()
{
    m_timers["generating linear systems"].start();
    generatePressureCorrectionSystem();
    m_timers["generating linear systems"].stop();

    m_timers["solving linear systems"].start();
    Field<Scalar> pCorrection = solveSystem(m_pressureSystemMatrix, m_pressureSystemSource);
    m_timers["solving linear systems"].stop();

    // Idk how, but explicit under relaxtion gives faster convergence than implicit
    pCorrection *= Config::pRelax;

    Field<Vector> uCorrection = getVelocityCorrection(pCorrection);
    Field<Scalar> massFluxesCorrection = getMassFluxesCorrection(pCorrection);

    m_pressureResidual = relativeResidual(m_currentPressure, pCorrection);

    m_currentVelocity          += uCorrection;
    m_currentPressure          += pCorrection;
    m_massFluxes += massFluxesCorrection;
}


void SimpleAlgorithm::initFields()
{
    Index totalCells = m_mesh.getCellAmount();
    Index totalFaces = m_mesh.getFaceAmount();

    // Initial guesses
    m_currentPressure           = Field<Scalar>::Constant(totalCells, 0);
    m_currentVelocity           = Field<Vector>::Constant(totalCells, {0,0,0});

    m_massFluxes = Field<Scalar>(totalFaces);
    m_pressureGradient      = Field<Vector>(totalCells);
    m_momentumSystemMatrix    = SparseMatrix(totalCells, totalCells);
    m_momentumSystemSource    = Matrix(totalCells, 3);
    m_pressureSystemMatrix    = SparseMatrix(totalCells, totalCells);
    m_pressureSystemSource    = Matrix(totalCells, 1);
    m_VbyA        = Field<Scalar>(totalCells);
    
    m_timers.clear();

    // Init mass fluxes
    // Can't be just zero because of boundary values
    auto uBoundaries = getVelocityBoundaries();
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        for (Index faceIdx : m_mesh.getCellFaces(cellIdx))
        {
            if
            (
                m_mesh.isBoundaryFace(faceIdx) 
                && m_mesh.getFaceBoundary(faceIdx).uBoundary.type == BoundaryConditionType::FIXED_VALUE
            )
            {
                m_massFluxes(faceIdx) = 
                (
                    Config::density *
                    m_mesh.getFaceBoundary(faceIdx).uBoundary.value.dot(m_mesh.getFaceVector(faceIdx))
                );
            }
            else
            {
                m_massFluxes(faceIdx) = 0;
            }
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
        m_pressureGradient(cellIdx) = 
        (
            Interpolation::computeCellGradient
            (
                m_mesh, cellIdx, pBoundaries, gradientScheme
            )
            .evaluate(m_currentPressure)
        );
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
        
        Vector faceVelocity = 
        (
            Interpolation::computeRhieChowVelocityOnFace
            (
                m_mesh, faceIdx, m_currentVelocity, m_currentPressure, m_pressureGradient, m_VbyA, uBoundaries, pBoundaries
            )
        );

        m_massFluxes(faceIdx) = faceVelocity.dot(faceVector) * Config::density;
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
        LinearCombination<Vector> convection = 
        (
            Interpolation::computeConvectionFluxOverCell
            (
                m_mesh, cellIdx, uBoundaries, m_massFluxes, convectionScheme
            )
        );
        
        LinearCombination<Vector> diffusion = 
        (
            Config::viscosity *
            Interpolation::computeDiffusionFluxOverCell
            (
                m_mesh, cellIdx, uBoundaries
            )
        );

        Vector pressureGradient = m_pressureGradient(cellIdx) * m_mesh.getCellVolume(cellIdx);

        LinearCombination<Vector> uEqn;
        uEqn += convection;
        uEqn -= diffusion;
        uEqn += pressureGradient;

        return uEqn;
    };

    generateSparseSystemImpl(m_momentumSystemMatrix, m_momentumSystemSource, m_mesh.getCellAmount(), uEqnGetter);
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

            Scalar VbyAf = 
            (
                Interpolation::valueOnFace
                (
                    m_mesh, faceIdx, zeroGradGetter<Scalar>()
                )
                .evaluate(m_VbyA)
            );
            
            diffusiveFlux +=
            (
                VbyAf * Config::density * faceVector.norm() *
                Interpolation::computeFaceNormalGradient
                (
                    m_mesh, cellIdx, faceIdx, pCorrBoundaries
                )
            );

            Scalar massFlux = m_massFluxes(faceIdx);
            if (cellIdx != m_mesh.getFaceOwner(faceIdx))
            {
                massFlux *= -1;
            }

            massFlow += massFlux;
        }

        LinearCombination<Scalar> pCorrEqn;
        pCorrEqn -= massFlow;
        pCorrEqn += diffusiveFlux;

        return pCorrEqn;
    };
    
    generateSparseSystemImpl(m_pressureSystemMatrix, m_pressureSystemSource, m_mesh.getCellAmount(), pCorrEqnGetter);
}


void SimpleAlgorithm::computeVbyA()
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto diagonal = m_momentumSystemMatrix.diagonal();

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        m_VbyA(cellIdx) = m_mesh.getCellVolume(cellIdx) / diagonal(cellIdx);
    }
    m_timers["explicit field computation"].stop();
}


Field<Vector> SimpleAlgorithm::getVelocityCorrection(Field<Scalar> const& pCorrection)
{
    m_timers["explicit field computation"].start();
    Index totalCells = m_mesh.getCellAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    Field<Vector> uCorrection(totalCells, 1);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        uCorrection(cellIdx) =
        (
            -m_VbyA(cellIdx) *
            Interpolation::computeCellGradient
            (
                m_mesh, cellIdx, pCorrBoundaries, gradientScheme
            )
            .evaluate(pCorrection)
        );
    }
    
    m_timers["explicit field computation"].stop();
    return uCorrection;
}


Field<Scalar> SimpleAlgorithm::getMassFluxesCorrection(Field<Scalar> const& pCorrection)
{
    m_timers["explicit field computation"].start();
    Index totalFaces = m_mesh.getFaceAmount();
    auto pCorrBoundaries = getPressureCorrectionBoundaries();
    Field<Scalar> massFluxesCorrection(totalFaces, 1);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Index ownerIdx = m_mesh.getFaceOwner(faceIdx);

        Scalar VbyAf =
        (
            Interpolation::valueOnFace
            (
                m_mesh, faceIdx, zeroGradGetter<Scalar>()
            )
            .evaluate(m_VbyA)
        );

        massFluxesCorrection(faceIdx) =
        (
            -Config::density * VbyAf *
            Interpolation::computeFaceNormalGradient
            (
                m_mesh, ownerIdx, faceIdx, pCorrBoundaries
            )
            .evaluate(pCorrection)
        );
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

    auto cmp = []<typename U>(Term<U> lhs, Term<U> rhs)
    {
        return lhs.idx < rhs.idx;
    };

#ifdef _OPENMP
    List<Triplet> threadTriplets;
    Index numThreads;
    #pragma omp parallel
    #pragma omp single
    {
        numThreads = omp_get_num_threads();
    }

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
            {
                threadTriplets.emplace_back(eqnIdx, varIdx, coeff);
            }

            rhs.row(eqnIdx) = -transpose(eqn.bias);
        }

        sizes[threadIdx] = threadTriplets.size();
        #pragma omp barrier
        
        #pragma omp single
        {
            for (Index idx = 1; idx <= numThreads; idx++)
            {
                prefix[idx] = prefix[idx-1] + sizes[idx-1];
            }
            
            triplets.resize(prefix.back());
        }

        for (Index idx = 0; idx < sizes[threadIdx]; idx++)
        {
            triplets[idx + prefix[threadIdx]] = threadTriplets[idx];
        }
    }
#else
    for (Index eqnIdx = 0; eqnIdx < size; eqnIdx++)
    {
        auto eqn = eqnGetter(eqnIdx);

        std::sort(eqn.terms.begin(), eqn.terms.end(), cmp);
        for (auto [coeff, varIdx] : eqn.terms)
        {
            triplets.emplace_back(eqnIdx, varIdx, coeff);
        }

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

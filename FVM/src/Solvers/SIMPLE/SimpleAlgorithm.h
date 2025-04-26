#pragma once

#include "Utils/Types.h"
#include "Utils/Timer.h"
#include "Mesh/MeshBase.h"
#include "Solvers/SolverBase.h"
#include "Discretization/Schemes/InterpolationSchemes.h"


class SimpleAlgorithm : public SolverBase
{
public:

    SimpleAlgorithm(MeshBase const& mesh);

    void solve() override;

    // Schemes
    Interpolation::Schemes::Gradient::Type gradientScheme;
    Interpolation::Schemes::Convection::Type convectionScheme;

private:

    // Fields in the current iteration
    Field<Vector> m_currentVelocity;
    Field<Scalar> m_currentPressure;
    Field<Vector> m_pressureGradient;
    Field<Scalar> m_VbyA;
    Field<Scalar> m_massFluxes;

    // Momentum matrix and rhs
    SparseMatrix m_momentumSystemMatrix;
    Matrix m_momentumSystemSource;

    // Pressure correction matrix and rhs
    SparseMatrix m_pressureSystemMatrix;
    Matrix m_pressureSystemSource;

    HashMap<std::string, Timer> m_timers;

    Scalar m_pressureResidual = 1;


    void initFields();
    void computePressureGradient();
    void solveMomentum();
    void computeMassFluxes();
    void correctPressure();

    bool converged() const;
    bool diverged() const;

    void generateMomentumSystem();
    void generatePressureCorrectionSystem();
    void computeVbyA();
    Field<Vector> getVelocityCorrection(Field<Scalar> const& pCorrection);
    Field<Scalar> getMassFluxesCorrection(Field<Scalar> const& pCorrection);

    Scalar relativeResidual(Matrix const& field, Matrix const& correction);

    BoundaryConditionGetter<Vector> getVelocityBoundaries();
    BoundaryConditionGetter<Scalar> getPressureBoundaries();
    BoundaryConditionGetter<Scalar> getPressureCorrectionBoundaries();
};

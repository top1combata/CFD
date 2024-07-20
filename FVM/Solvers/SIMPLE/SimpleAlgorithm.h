#pragma once

#include "Utils/Types.h"
#include "Utils/Timer.h"
#include "Mesh/MeshBase.h"


class SimpleAlgorithm
{
public:

    SimpleAlgorithm(MeshBase const& mesh);

    void solve();

    VectorField& getU();
    ScalarField& getP();

private:

    // fields
    VectorField m_U;
    ScalarField m_p;
    VectorField m_p_grad;
    ScalarField m_VbyA;
    ScalarField m_mass_fluxes;
    // momentum matrix and
    SparseMatrix m_U_matrix;
    Matrix m_U_source;
    // pressure correction matrix
    SparseMatrix m_p_matrix;
    Matrix m_p_source;

    HashMap<std::string, Timer> m_timers;

    Scalar m_U_residual = 1;
    Scalar m_p_residual = 1;

    // associated mesh
    MeshBase const& m_mesh;

    void initFields();
    void computePressureGradient();
    void solveMomentum();
    void computeMassFluxes();
    void correctPressure();

    bool converged();

    void generateMomentumSystem();
    void generatePressureCorrectionSystem();
    void computeVbyA();
    VectorField getVelocityCorrection(ScalarField const& pCorrection);
    ScalarField getMassFluxesCorrection(ScalarField const& pCorrection);


    Scalar relativeResidual(Matrix const& field, Matrix const& correction);

    BoundaryConditionGetter<Vector> getVelocityBoundaries();
    BoundaryConditionGetter<Scalar> getPressureBoundaries();
    BoundaryConditionGetter<Scalar> getPressureCorrectionBoundaries();
};
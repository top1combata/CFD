#pragma once

#include "Utils/Types.h"
#include "Utils/Timer.h"
#include "Mesh/MeshBase.h"


class SimpleAlgorithm
{
public:

    SimpleAlgorithm(MeshBase const& mesh);

    void solve();

    Field<Vector>& getU();
    Field<Scalar>& getP();

private:

    // fields
    Field<Vector> m_U;
    Field<Scalar> m_p;
    Field<Vector> m_p_grad;
    Field<Scalar> m_VbyA;
    Field<Scalar> m_mass_fluxes;
    // momentum matrix and
    SparseMatrix m_U_matrix;
    Matrix m_U_source;
    // pressure correction matrix
    SparseMatrix m_p_matrix;
    Matrix m_p_source;

    HashMap<std::string, Timer> m_timers;

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
    Field<Vector> getVelocityCorrection(Field<Scalar> const& pCorrection);
    Field<Scalar> getMassFluxesCorrection(Field<Scalar> const& pCorrection);


    Scalar relativeResidual(Matrix const& field, Matrix const& correction);

    BoundaryConditionGetter<Vector> getVelocityBoundaries();
    BoundaryConditionGetter<Scalar> getPressureBoundaries();
    BoundaryConditionGetter<Scalar> getPressureCorrectionBoundaries();
};
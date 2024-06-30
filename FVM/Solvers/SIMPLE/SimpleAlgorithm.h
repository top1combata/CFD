#pragma once

#include "Utils/Types.h"
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
    // momentum matrix
    Matrix m_Au;

    Scalar m_U_residual = 1;
    Scalar m_p_residual = 1;

    // associated mesh
    MeshBase const&m_mesh;

    void initFields();
    void computePressureGradient();
    void solveMomentum();
    void updateMassFluxes();
    void correctPressure();
    bool converged();

    static Scalar relativeResidual(Matrix const& field, Matrix const& correction);

    BoundaryConditionGetter<Vector> uBoundaries();
    BoundaryConditionGetter<Scalar> pBoundaries();
};
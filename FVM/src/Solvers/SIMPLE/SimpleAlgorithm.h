#pragma once

#include "Utils/Types.h"
#include "Utils/Timer.h"
#include "Mesh/MeshBase.h"
#include "Discretization/Schemes/InterpolationSchemes.h"


class SimpleAlgorithm
{
public:

    SimpleAlgorithm(MeshBase const& mesh);

    void solve();

    Field<Vector>& getU();
    Field<Scalar>& getP();

    // Schemes
    Interpolation::Schemes::Gradient::Type gradientScheme;
    Interpolation::Schemes::Convection::Type convectionScheme;

private:

    // fields
    Field<Vector> m_U;
    Field<Scalar> m_p;
    Field<Vector> m_pGrad;
    Field<Scalar> m_VbyA;
    Field<Scalar> m_massFluxes;
    // momentum matrix and
    SparseMatrix m_UMatrix;
    Matrix m_USource;
    // pressure correction matrix
    SparseMatrix m_pMatrix;
    Matrix m_pSource;

    HashMap<std::string, Timer> m_timers;

    Scalar m_pResidual = 1;

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

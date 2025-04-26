#include "SolverBase.h"


SolverBase::SolverBase(MeshBase const& mesh) : m_mesh(mesh) {}

MeshBase const& SolverBase::getMesh() const
{
    return m_mesh;
}

Field<Vector> const& SolverBase::getVelocity(Index timePointIdx) const
{
    return m_velocity[timePointIdx];
}

Field<Scalar> const& SolverBase::getPressure(Index timePointIdx) const
{
    return m_pressure[timePointIdx];
}

bool SolverBase::isConverged() const
{
    return m_converged;
}

#include "SolverBase.h"

#include <cassert>


SolverBase::SolverBase(MeshBase const& mesh) : m_mesh(mesh) {}

bool SolverBase::isTransient() const
{
    return m_isTransient;
}

void SolverBase::setTransient(bool state)
{
    m_isTransient = state;
}

Index SolverBase::getTimePointAmount() const
{
    assert(m_velocity.size() == m_pressure.size());
    return m_velocity.size();
}

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

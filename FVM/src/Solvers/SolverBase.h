#pragma once

#include "Utils/Types.h"
#include "Mesh/MeshBase.h"


class SolverBase
{
public:

    SolverBase(MeshBase const& mesh);

    virtual void solve() = 0;
    virtual bool isConverged() const;

    MeshBase const& getMesh() const;
    Field<Vector> const& getVelocity(Index timePointIdx) const;
    Field<Scalar> const& getPressure(Index timePointIdx) const;

protected:

    MeshBase const& m_mesh;
    List<Field<Vector>> m_velocity;
    List<Field<Scalar>> m_pressure;

    bool m_converged;
};

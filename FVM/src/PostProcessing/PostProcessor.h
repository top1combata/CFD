#pragma once

#include "Solvers/SolverBase.h"


class PostProcessor
{
public:

    PostProcessor(SolverBase const& solver);

    void show();

private:

    MeshBase const& m_mesh;
    SolverBase const& m_solver;
};

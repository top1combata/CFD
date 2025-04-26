#pragma once

#include "Primitives.h"
#include "Mesh/MeshBase.h"
#include "Solvers/SolverBase.h"


sf::View getInitializeView(MeshBase const& mesh);

List<sf::ConvexShape> getMeshCells(MeshBase const& mesh);

List<Line> getMeshFaces(MeshBase const& mesh);

List<Arrow> getVelocityField(SolverBase const& solver);

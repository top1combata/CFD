#pragma once

#include "Primitives.h"
#include "Mesh/MeshBase.h"
#include "Solvers/SolverBase.h"


sf::View getInitializeView(MeshBase const& mesh);

List<sf::ConvexShape> getMeshCells(MeshBase const& mesh);

List<Line> getMeshFaces(MeshBase const& mesh);

List<Arrow> getVelocityField(SolverBase const& solver, Index timePointIdx, bool useAbsoluteVectorLength);

// Argument should be in [0, 1]
sf::Color getGradientColor(Scalar t);

void heatMapColor(List<sf::ConvexShape>& cells, std::function<Scalar(Index)> const& fieldGetter);

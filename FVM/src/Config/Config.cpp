#include "Config.h"

#include <limits>


Scalar Config::density = 1e3;
Scalar Config::viscosity = 1e-3;

Scalar Config::uRelax = 0.7;
Scalar Config::pRelax = 0.3;

Scalar Config::uTolerance = 1e-5;
Scalar Config::pTolerance = 1e-5;
Index  Config::maxIterations = 10'000;

Scalar Config::uSystemTolerance = std::numeric_limits<Scalar>::epsilon();
Scalar Config::pSystemTolerance = std::numeric_limits<Scalar>::epsilon();

Scalar Config::timeStep = 0;
Scalar Config::timeBegin = 0;
Scalar Config::timeEnd = 0;

using namespace Interpolation::Schemes;
Gradient::Type Config::gradientScheme = Gradient::GREEN_GAUSE;
Convection::Type Config::convectionScheme = Convection::SOU;

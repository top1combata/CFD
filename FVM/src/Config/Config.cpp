#include "Config.h"

Scalar Config::density = 1e3;
Scalar Config::viscosity = 1e-3;

Scalar Config::uRelax = 0.7;
Scalar Config::pRelax = 0.3;

Scalar Config::uTolerance = 1e-5;
Scalar Config::pTolerance = 1e-5;
Index  Config::maxIterations = 10'000;
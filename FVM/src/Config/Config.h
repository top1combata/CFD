#include "Utils/Types.h"
#include "Discretization/Schemes/InterpolationSchemes.h"


class Config
{
public:

    static Scalar density;
    // Dynamic viscosity of fluid
    static Scalar viscosity;
    
    static Scalar uRelax;
    static Scalar pRelax;

    static Scalar uTolerance;
    static Scalar pTolerance;
    static Index maxIterations;

    // Interpolation schemes
    static Interpolation::Schemes::Gradient::Type gradientScheme;
    static Interpolation::Schemes::Convection::Type convectionScheme;
};

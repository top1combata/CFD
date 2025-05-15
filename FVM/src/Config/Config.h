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

    // For transient solver
    static Scalar timeStep;
    static Scalar timeBegin;
    static Scalar timeEnd;

    static Scalar uTolerance;
    static Scalar pTolerance;
    static Index maxIterations;
    
    // Tolerance for lienar system solvers
    static Scalar uSystemTolerance;
    static Scalar pSystemTolerance;

    // Interpolation schemes
    static Interpolation::Schemes::Gradient::Type gradientScheme;
    static Interpolation::Schemes::Convection::Type convectionScheme;
};

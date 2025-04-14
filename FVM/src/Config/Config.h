#include "Utils/Types.h"


class Config
{
public:

    static Scalar density;
    // dynamic viscosity of fluid
    static Scalar viscosity;
    
    static Scalar uRelax;
    static Scalar pRelax;

    static Scalar uTolerance;
    static Scalar pTolerance;
    static Index maxIterations;
};

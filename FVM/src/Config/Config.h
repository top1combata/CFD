#include "Utils/Types.h"


class Config
{
public:

    static Scalar density;
    static Scalar viscosity;
    
    static Scalar uRelax;
    static Scalar pRelax;

    static Scalar uTolerance;
    static Scalar pTolerance;
    static Index maxIterations;
};
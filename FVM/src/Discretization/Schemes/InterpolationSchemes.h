#pragma once


namespace Interpolation::Schemes::Gradient
{

enum Type
{
    GREEN_GAUSE,
    LEAST_SQAURE
};

}

namespace Interpolation::Schemes::Convection
{

enum Type
{
    CENTRAL_DIFFERENCE,
    UPWIND,
    DOWNWIND,
    FROMM,
    SOU,
    QUICK
};

}

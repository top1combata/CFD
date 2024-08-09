#pragma once

#include "Types.h"


inline Scalar getFieldValue(ScalarField const& field, Index idx) {return field(idx);}

inline Vector getFieldValue(VectorField const& field, Index idx) {return field.row(idx).transpose();}

template<class T>
T zero() {return static_cast<T>(0);}

template<class T>
T zero()
requires requires {{T::Zero()} -> std::convertible_to<T>;}
{return T::Zero();}
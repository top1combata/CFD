#pragma once

#include "Types.h"


template<class T>
T zero() {return static_cast<T>(0);}

template<class T>
T zero()
requires requires {{T::Zero()} -> std::convertible_to<T>;}
{return T::Zero();}


template<class U, class V>
struct ProductType;

template<class U>
struct ProductType<U, Scalar> {using type = U;};

template<class U>
struct ProductType<Scalar, U> {using type = U;};

template<>
struct ProductType<Scalar, Scalar> {using type = Scalar;};

template<>
struct ProductType<Vector, Vector> {using type = Tensor;};
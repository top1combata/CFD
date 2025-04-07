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


inline Scalar innerProduct(Vector lhs, Vector rhs)
{
    return lhs.dot(rhs);
}

inline Vector innerProduct(Tensor const& lhs, Vector rhs)
{
    return lhs * rhs;
}

inline Vector innerProduct(Vector lhs, Tensor const& rhs)
{
    return (lhs.transpose() * rhs).transpose();
}

inline Tensor innerProduct(Tensor const& lhs, Tensor const& rhs)
{
    return lhs * rhs;
}
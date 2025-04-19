#pragma once

#include "Utils/Types.h"
#include "Utils/TypesOperations.h"

template<class CoeffType>
class Term
{
public:

    CoeffType coeff;
    Index idx;
};


template<class VarType, class CoeffType = Scalar>
class LinearCombination
{
public:

    using BiasType = ProductType<VarType, CoeffType>::type;

    List<Term<CoeffType>> terms;
    BiasType bias = zero<BiasType>();

    LinearCombination() = default;

    explicit LinearCombination(BiasType const& value) : bias(value) {}

    LinearCombination(std::initializer_list<Term<CoeffType>> const&);

    LinearCombination& operator+=(BiasType const& value);

    LinearCombination& operator-=(BiasType const& value);

    LinearCombination& operator+=(Term<CoeffType> const&);

    LinearCombination& operator-=(Term<CoeffType> const&);

    LinearCombination& operator+=(LinearCombination const&);

    LinearCombination& operator-=(LinearCombination const&);

    LinearCombination& operator*=(Scalar);

    LinearCombination& operator/=(Scalar);

    BiasType evaluate(Field<VarType> const&) const;

    LinearCombination<VarType, Scalar> dot(Vector) const
    requires std::same_as<CoeffType, Vector>;
};

template<class U, class V>
LinearCombination<U,V> operator+(LinearCombination<U,V>, LinearCombination<U,V> const&);

template<class U, class V>
LinearCombination<U,V> operator+(LinearCombination<U,V>, typename LinearCombination<U,V>::BiasType const&);

template<class U, class V>
LinearCombination<U,V> operator+(typename LinearCombination<U,V>::BiasType const&, LinearCombination<U,V>);

template<class U, class V>
LinearCombination<U,V> operator-(LinearCombination<U,V>, LinearCombination<U,V> const&);

template<class U, class V>
LinearCombination<U,V> operator-(LinearCombination<U,V>, typename LinearCombination<U,V>::BiasType const&);

template<class U, class V>
LinearCombination<U,V> operator-(typename LinearCombination<U,V>::BiasType const&, LinearCombination<U,V>);

template<class U, class V>
LinearCombination<U,V> operator-(LinearCombination<U,V>);

template<class U, class V>
LinearCombination<U,V> operator/(LinearCombination<U,V>, Scalar);


// Product with scalar or vector
// In case of multiplying by vector coeffs of the result calculated as outer prodcut
template<class VarType, class CoeffType, class T>
LinearCombination<VarType, typename ProductType<CoeffType, T>::type>
operator*(LinearCombination<VarType, CoeffType> const&, T)
requires std::same_as<T, Scalar> || std::same_as<T, Vector>;

template<class VarType, class CoeffType, class T>
LinearCombination<VarType, typename ProductType<T, CoeffType>::type>
operator*(T, LinearCombination<VarType, CoeffType> const&)
requires std::same_as<T, Scalar> || std::same_as<T, Vector>;

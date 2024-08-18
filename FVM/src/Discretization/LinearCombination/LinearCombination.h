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
LinearCombination<U,V> operator*(LinearCombination<U,V>, Scalar);

template<class U, class V>
LinearCombination<U,V> operator*(Scalar, LinearCombination<U,V>);

template<class U, class V>
LinearCombination<U,V> operator/(LinearCombination<U,V>, Scalar);


LinearCombination<Scalar, Vector> operator*(LinearCombination<Scalar, Scalar> const&, Vector);
LinearCombination<Scalar, Vector> operator*(Vector, LinearCombination<Scalar, Scalar> const&);

LinearCombination<Vector, Vector> operator*(LinearCombination<Vector, Scalar> const&, Vector);
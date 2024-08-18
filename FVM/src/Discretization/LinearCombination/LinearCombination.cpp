#include "LinearCombination.h"
#include <algorithm>


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>::LinearCombination(std::initializer_list<Term<CoeffType>> const& list) : terms(list.size())
{
    std::copy(list.begin(), list.end(), terms.begin());
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator+=(BiasType const& value)
{
    bias += value;
    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator-=(BiasType const& value)
{
    bias -= value;
    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator+=(Term<CoeffType> const& term)
{
    auto iter = std::find_if(terms.begin(), terms.end(),
    [&term](Term<CoeffType> const& t) 
    {
        return t.idx == term.idx;
    });
    if (iter == terms.end())
        terms.push_back(term);
    else
        iter->coeff += term.coeff;
    
    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator-=(Term<CoeffType> const& term)
{
    Term<CoeffType> copy = term;
    copy.coeff *= -1;
    *this += copy;
    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator+=(LinearCombination const& rhs)
{
    bias += rhs.bias;
    for (auto const& term : rhs.terms)
        *this += term;

    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator-=(LinearCombination const& rhs)
{
    bias -= rhs.bias;
    for (auto const& term : rhs.terms)
        *this -= term;

    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator*=(Scalar scalar)
{
    bias *= scalar;
    for (auto& term : terms)
        term.coeff *= scalar;
    
    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator/=(Scalar scalar)
{
    return *this *= 1/scalar;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>::BiasType LinearCombination<VarType, CoeffType>::evaluate(Field<VarType> const& field) const
{
    BiasType res = bias;
    for (auto const& [coeff, idx] : terms)
        res += field(idx) * coeff;

    return res;
}


template<>
Tensor LinearCombination<Vector, Vector>::evaluate(Field<Vector> const& field) const
{
    Tensor res = bias;
    for (auto [coeff, idx] : terms)
        res += field(idx) * coeff.transpose();

    return res;
}


template<>
LinearCombination<Scalar, Scalar> LinearCombination<Scalar, Vector>::dot(Vector vec) const
{
    LinearCombination<Scalar, Scalar> res(bias.dot(vec));
    res.terms.reserve(terms.size());

    for (auto [coeff, idx] : terms)
        res.terms.emplace_back(coeff.dot(vec), idx);
    
    return res;
}


template<>
LinearCombination<Vector, Scalar> LinearCombination<Vector, Vector>::dot(Vector vec) const
{
    LinearCombination<Vector, Scalar> res(bias * vec);
    res.terms.reserve(terms.size());

    for (auto [coeff, idx] : terms)
        res.terms.emplace_back(coeff.dot(vec), idx);
    
    return res;
}


template<class U, class V>
LinearCombination<U,V> operator+(LinearCombination<U,V> lhs, LinearCombination<U,V> const& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U,V> operator+(LinearCombination<U,V> lhs, typename LinearCombination<U,V>::BiasType const& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U,V> operator+(typename LinearCombination<U,V>::BiasType const& lhs, LinearCombination<U,V> rhs)
{
    rhs += lhs;
    return rhs;
}

template<class U, class V>
LinearCombination<U,V> operator-(LinearCombination<U,V> lhs, LinearCombination<U,V> const& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U,V> operator-(LinearCombination<U,V> lhs, typename LinearCombination<U,V>::BiasType const& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U,V> operator-(typename LinearCombination<U,V>::BiasType const& lhs, LinearCombination<U,V> rhs)
{
    rhs -= lhs;
    rhs *= -1;
    return rhs;
}

template<class U, class V>
LinearCombination<U,V> operator-(LinearCombination<U,V> lc)
{
    lc *= -1;
    return lc;
}

template<class U, class V>
LinearCombination<U,V> operator*(LinearCombination<U,V> lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U,V> operator*(Scalar lhs, LinearCombination<U,V> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<class U, class V>
LinearCombination<U,V> operator/(LinearCombination<U,V> lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}


LinearCombination<Scalar, Vector> operator*(LinearCombination<Scalar, Scalar> const& lc, Vector vec)
{
    LinearCombination<Scalar, Vector> res(vec * lc.bias);
    res.terms.reserve(lc.terms.size());
    for (auto [coeff, idx] : lc.terms)
        res.terms.emplace_back(coeff*vec, idx);
    return res;
}


LinearCombination<Scalar, Vector> operator*(Vector vec, LinearCombination<Scalar, Scalar> const& lc)
{
    auto res = lc*vec;
    return res;
}


LinearCombination<Vector, Vector> operator*(LinearCombination<Vector, Scalar> const& lc, Vector vec)
{
    LinearCombination<Vector, Vector> res(lc.bias * vec.transpose());
    res.terms.reserve(lc.terms.size());
    for (auto [coeff, idx] : lc.terms)
        res.terms.emplace_back(coeff*vec, idx);
    return res;
}


// explicit instantiation non class members operators
// using a macro to avoid repeating
#define INSTANTIATE(U,V)                                                                                            \
template class LinearCombination<U, V>;                                                                             \
template LinearCombination<U,V> operator+(LinearCombination<U,V>, LinearCombination<U,V> const&);                   \
template LinearCombination<U,V> operator+(LinearCombination<U,V>, typename LinearCombination<U,V>::BiasType const&);\
template LinearCombination<U,V> operator+(typename LinearCombination<U,V>::BiasType const&, LinearCombination<U,V>);\
template LinearCombination<U,V> operator-(LinearCombination<U,V>, LinearCombination<U,V> const&);                   \
template LinearCombination<U,V> operator-(LinearCombination<U,V>, typename LinearCombination<U,V>::BiasType const&);\
template LinearCombination<U,V> operator-(typename LinearCombination<U,V>::BiasType const&, LinearCombination<U,V>);\
template LinearCombination<U,V> operator-(LinearCombination<U,V>);                                                  \
template LinearCombination<U,V> operator*(LinearCombination<U,V>, Scalar);                                          \
template LinearCombination<U,V> operator*(Scalar, LinearCombination<U,V>);                                          \
template LinearCombination<U,V> operator/(LinearCombination<U,V>, Scalar);

INSTANTIATE(Scalar, Scalar)
INSTANTIATE(Scalar, Vector)
INSTANTIATE(Vector, Scalar)
INSTANTIATE(Vector, Vector)
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
    {
        terms.push_back(term);
    }
    else
    {
        iter->coeff += term.coeff;
    }
    
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
    {
        *this += term;
    }

    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator-=(LinearCombination const& rhs)
{
    bias -= rhs.bias;
    for (auto const& term : rhs.terms)
    {
        *this -= term;
    }

    return *this;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, CoeffType>& LinearCombination<VarType, CoeffType>::operator*=(Scalar scalar)
{
    bias *= scalar;
    for (auto& term : terms)
    {
        term.coeff *= scalar;
    }
    
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
    for (auto [coeff, idx] : terms)
    {
        res += outerProduct(field(idx), coeff);
    }

    return res;
}


template<class VarType, class CoeffType>
LinearCombination<VarType, Scalar> LinearCombination<VarType, CoeffType>::dot(Vector vec) const
requires std::same_as<CoeffType, Vector>
{
    LinearCombination<VarType, Scalar> res(innerProduct(bias, vec));
    res.terms.reserve(terms.size());

    for (auto [coeff, idx] : terms)
    {
        res.terms.emplace_back(coeff.dot(vec), idx);
    }
    
    return res;
}


template<class U, class V>
LinearCombination<U, V> operator+(LinearCombination<U, V> lhs, LinearCombination<U, V> const& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U, V> operator+(LinearCombination<U, V> lhs, typename LinearCombination<U, V>::BiasType const& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U, V> operator+(typename LinearCombination<U, V>::BiasType const& lhs, LinearCombination<U, V> rhs)
{
    rhs += lhs;
    return rhs;
}

template<class U, class V>
LinearCombination<U, V> operator-(LinearCombination<U, V> lhs, LinearCombination<U, V> const& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U, V> operator-(LinearCombination<U, V> lhs, typename LinearCombination<U, V>::BiasType const& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U, V> operator-(typename LinearCombination<U, V>::BiasType const& lhs, LinearCombination<U, V> rhs)
{
    rhs -= lhs;
    rhs *= -1;
    return rhs;
}

template<class U, class V>
LinearCombination<U, V> operator-(LinearCombination<U, V> lc)
{
    lc *= -1;
    return lc;
}

template<class U, class V>
LinearCombination<U, V> operator*(LinearCombination<U, V> lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class U, class V>
LinearCombination<U, V> operator*(Scalar lhs, LinearCombination<U, V> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<class U, class V>
LinearCombination<U, V> operator/(LinearCombination<U, V> lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}


template<class VarType, class CoeffType, class T>
LinearCombination<VarType, typename ProductType<CoeffType, T>::type>
operator*(LinearCombination<VarType, CoeffType> const& lc, T value)
requires std::same_as<T, Scalar> || std::same_as<T, Vector>
{
    LinearCombination<VarType, typename ProductType<CoeffType, T>::type> res;
    
    res.bias = (outerProduct(lc.bias, value));
    res.terms.reserve(lc.terms.size());
    for (auto [coeff, idx] : lc.terms)
    {
        res.terms.emplace_back(outerProduct(coeff, value), idx);
    }

    return res;
}


template<class VarType, class CoeffType, class T>
LinearCombination<VarType, typename ProductType<T, CoeffType>::type>
operator*(T value, LinearCombination<VarType, CoeffType> const& lc)
requires std::same_as<T, Scalar> || std::same_as<T, Vector>
{
    LinearCombination<VarType, typename ProductType<T, CoeffType>::type> res;

    // implementing not using previous, cause of order of outer product
    res.bias = (outerProduct(value, lc.bias));
    res.terms.reserve(lc.terms.size());
    for (auto [coeff, idx] : lc.terms)
    {
        res.terms.emplace_back(outerProduct(value, coeff), idx);
    }

    return res;
}


// explicit instantiation non class members operators
// using a macro to avoid repeating
#define INSTANTIATE_COMMON(U, V)                                                                                        \
template class LinearCombination<U, V>;                                                                                 \
template LinearCombination<U, V> operator+(LinearCombination<U, V>, LinearCombination<U, V> const&);                    \
template LinearCombination<U, V> operator+(LinearCombination<U, V>, typename LinearCombination<U, V>::BiasType const&); \
template LinearCombination<U, V> operator+(typename LinearCombination<U, V>::BiasType const&, LinearCombination<U, V>); \
template LinearCombination<U, V> operator-(LinearCombination<U, V>, LinearCombination<U, V> const&);                    \
template LinearCombination<U, V> operator-(LinearCombination<U, V>, typename LinearCombination<U, V>::BiasType const&); \
template LinearCombination<U, V> operator-(typename LinearCombination<U, V>::BiasType const&, LinearCombination<U, V>); \
template LinearCombination<U, V> operator-(LinearCombination<U, V>);                                                    \
template LinearCombination<U, V> operator/(LinearCombination<U, V>, Scalar);                                            \
template LinearCombination<U, ProductType<V, Scalar>::type> operator*(LinearCombination<U, V> const&, Scalar);          \
template LinearCombination<U, ProductType<Scalar, V>::type> operator*(Scalar, LinearCombination<U, V> const&);          \

#define INSTANTIATE_MULTIPLY_BY_VECTOR(U, V)                                                                            \
template LinearCombination<U, ProductType<V, Vector>::type> operator*(LinearCombination<U, V> const&, Vector);          \
template LinearCombination<U, ProductType<Vector, V>::type> operator*(Vector, LinearCombination<U, V> const&);          \

INSTANTIATE_COMMON(Scalar, Scalar)
INSTANTIATE_COMMON(Scalar, Vector)
INSTANTIATE_COMMON(Vector, Scalar)
INSTANTIATE_COMMON(Vector, Vector)

INSTANTIATE_MULTIPLY_BY_VECTOR(Scalar, Scalar)
INSTANTIATE_MULTIPLY_BY_VECTOR(Scalar, Vector)
INSTANTIATE_MULTIPLY_BY_VECTOR(Vector, Scalar)
// not instantiating for <Vector, Vector> cause we don't use tensors of 3rd rank or higher

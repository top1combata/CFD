#pragma once

#include "Utils/Types.h"

class Term
{
public:

    Scalar coeff;
    Index idx;

    bool operator==(Term another) const
    {
        return idx == another.idx;
    }
};

template<class T>
class LinearCombination
{
public:

    List<Term> terms;
    T bias;

    LinearCombination(T value);

    LinearCombination(std::initializer_list<Term> const& lst);


    LinearCombination& operator+=(LinearCombination const& rhs);

    LinearCombination& operator-=(LinearCombination const& rhs);
};


template<class T>
LinearCombination<T> operator+(LinearCombination<T> lhs, LinearCombination<T> const& rhs);

template<class T>
LinearCombination<T> operator-(LinearCombination<T> lhs, LinearCombination<T> const& rhs);

template<class T>
LinearCombination<T> operator-(LinearCombination<T> lc);

template<class T>
LinearCombination<T>& operator*=(LinearCombination<T>& lhs, Scalar rhs);

template<class T>
LinearCombination<T> operator*(LinearCombination<T> lhs, Scalar rhs);

template<class T>
LinearCombination<T> operator*(Scalar lhs, LinearCombination<T> rhs);

template<class T>
std::ostream& operator<<(std::ostream& os, LinearCombination<T> const& lc);

#include "LinearCombination.hpp"
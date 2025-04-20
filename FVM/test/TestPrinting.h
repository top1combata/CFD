#pragma once

#include <Utils/Types.h>
#include <Discretization/LinearCombination.h>


inline void PrintTo(Scalar scalar, std::ostream* os)
{
    *os << scalar;
}


inline void PrintTo(Vector vec, std::ostream* os)
{
    *os << "Vector(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
}


inline void PrintTo(Tensor const& tensor, std::ostream* os)
{
    *os << "Tensor(";
    for (Index i = 0; i < 3; i++)
    {
        *os << "{" << tensor(i, 0) << ", " << tensor(i, 1) << ", " << tensor(i, 2) << "}";
        if (i < 2)
        {
            *os << ", ";
        }
    }
}


template<class U>
void PrintTo(Term<U> const& term, std::ostream* os)
{
    *os << "Term{Coeff(";
    PrintTo(term.coeff, os);
    *os << "), " << "Idx(" << term.idx << ")}";
}


template<class U, class V>
void PrintTo(LinearCombination<U, V> const& lc, std::ostream* os)
{
    *os << "LinearCombination{" << "Bias(";
    PrintTo(lc.bias, os);
    *os << "), ";
    for (auto const& term : lc.terms)
    {
        PrintTo(term, os);
        *os << ", ";
    }
    *os << "}";
}

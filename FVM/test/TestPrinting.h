#pragma once

#include <Utils/Types.h>
#include <Boundary/BoundaryCondition.h>
#include <Discretization/LinearCombination.h>

#include <sstream>


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


using std::to_string;

inline std::string to_string(Vector vec)
{
    std::stringstream ss;
    PrintTo(vec, &ss);
    return ss.str();
}


inline std::string to_string(Tensor tensor)
{
    std::stringstream ss;
    PrintTo(tensor, &ss);
    return ss.str();
}


template<class T>
std::string to_string(Term<T> const& term)
{
    std::stringstream ss;
    PrintTo(term, &ss);
    return ss.str();
}


template<class U, class V>
std::string to_string(LinearCombination<U, V> const& lc)
{
    std::stringstream ss;
    PrintTo(lc, &ss);
    return ss.str();
}


inline std::string to_string(BoundaryConditionType type)
{
    return (type == BoundaryConditionType::FIXED_VALUE ? "fixed value" : "fixed gradient");
}


template <class T> std::string to_string(BoundaryCondition<T> const& boundary)
{
    return std::string("{")
        .append(to_string(boundary.type))
        .append(", ")
        .append(to_string(boundary.value))
        .append("}")
    ;
}


inline std::string to_string(Boundaries const& boundaries)
{
    return std::string("{U")
        .append(to_string(boundaries.uBoundary))
        .append(", p")
        .append(to_string(boundaries.pBoundary))
        .append("}")
    ;
}

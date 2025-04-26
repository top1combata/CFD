#pragma once

#include "TestPrinting.h"
#include <gtest/gtest.h>

#include <Boundary/BoundaryCondition.h>
#include <Discretization/LinearCombination.h>
#include <Utils/Types.h>

#include <random>


#define EXPECT_MAP_VALUE(map, key, value)                                   \
{                                                                           \
    EXPECT_TRUE((map).contains(key) && (map)[key] == (value))               \
        << ((map).contains(key)                                             \
        ? std::string("value is not exect at key ").append(to_string(key))  \
            .append(", returned ").append(to_string((map)[key]))            \
            .append(", expected ").append(to_string(value))                 \
        : std::string("no value at key ").append(to_string(key)));          \
}

inline bool operator==(BoundaryCondition<Scalar> const& b1, BoundaryCondition<Scalar> const& b2)
{
    return b1.type == b2.type && abs(b1.value - b2.value) / std::max(abs(b1.value), Scalar(1)) < 1e-7;
}

inline bool operator==(BoundaryCondition<Vector> const& b1, BoundaryCondition<Vector> const& b2)
{
    return b1.type == b2.type && (b1.value - b2.value).norm() / std::max(b1.value.norm(), Scalar(1)) < 1e-7;
}

inline bool operator==(Boundaries const& b1, Boundaries const& b2)
{
    return b1.uBoundary == b2.uBoundary && b1.pBoundary == b2.pBoundary;
}

template <class T>
    requires std::convertible_to<T, Scalar>
Scalar norm(T value)
{
    return std::abs(value);
}

template <class T>
    requires requires(T value) {
        { value.norm() };
    }
Scalar norm(T const& value)
{
    return value.norm();
}

template <class T> bool operator==(Term<T> const& t1, Term<T> const& t2)
{
    constexpr Scalar tolerance = 1e-5;
    return t1.idx == t2.idx && norm(t1.coeff - t2.coeff) < tolerance;
}

template <class U, class V> bool operator==(LinearCombination<U, V> const& lc1, LinearCombination<U, V> const& lc2)
{
    constexpr Scalar tolerance = 1e-5;
    auto             lc = lc1 - lc2;
    bool             res = norm(lc.bias) < tolerance;
    for (auto const& [coeff, _] : lc.terms)
        if (norm(coeff) > tolerance)
            res = false;
    return res;
}

template <class T> std::ostream& operator<<(std::ostream& os, Term<T> const& term)
{
    return os << "{" << term.coeff << "," << term.idx << "}";
}

template <class U, class V> void sortIndices(LinearCombination<U, V>& lc)
{
    std::sort(lc.terms.begin(), lc.terms.end(), [](auto t1, auto t2) { return t1.idx < t2.idx; });
}

// Generation of entities

inline Index randomIndex(Index min, Index max)
{
    constexpr uint64_t SEED = 42;

    static std::mt19937_64 gen(SEED);
    int64_t                range = max - min + 1;

    return gen() % range + min;
}

inline Scalar randomScalar()
{
    Index num = randomIndex(-100, 100);
    Index denum = (1 << randomIndex(0, 5));
    return num / (Scalar)denum;
}

inline Vector randomVector()
{
    Vector res;
    for (int idx = 0; idx < 3; idx++)
        res(idx) = randomScalar();
    return res;
}

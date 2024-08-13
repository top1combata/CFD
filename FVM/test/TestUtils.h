#pragma once

#include <gtest/gtest.h>
#include <format>
#include <Utils/Types.h>
#include <Boundary/BoundaryCondition.h>


using std::to_string;

static std::string to_string(Vector const& vec)
{
    return std::format("({}, {}, {})", vec.x(), vec.y(), vec.z());
}

static std::string to_string(BoundaryConditionType type)
{
    return (type == BoundaryConditionType::FIXED_VALUE ? "fixed value" : "fixed gradient");
}

template<class T>
std::string to_string(BoundaryCondition<T> const& boundary)
{
    return std::format("{{{} {}}}", to_string(boundary.type), to_string(boundary.value));
}

static std::string to_string(Boundaries const& boundaries)
{
    return std::format("{{U{} , p{}}}", to_string(boundaries.uBoundary), to_string(boundaries.pBoundary));
}

#define EXPECT_MAP_VALUE(map, key, value) {                                     \
    EXPECT_TRUE((map).contains(key) && (map)[key] == (value)) <<                \
        ((map).contains(key) ?                                                  \
        std::format("value is not exected at key {}, returned {}, expected {}", \
            key, to_string((map)[key]), to_string(value)) :                     \
        std::format("no value at key {}", key));                                \
}


static bool operator==(BoundaryCondition<Scalar> const& b1, BoundaryCondition<Scalar> const& b2)
{
    return b1.type == b2.type && abs(b1.value - b2.value) / std::max(abs(b1.value), Scalar(1)) < 1e-7;
}

static bool operator==(BoundaryCondition<Vector> const& b1, BoundaryCondition<Vector> const& b2)
{
    return b1.type == b2.type && (b1.value - b2.value).norm() / std::max(b1.value.norm(), Scalar(1)) < 1e-7;
}

static bool operator==(Boundaries const& b1, Boundaries const& b2)
{
    return b1.uBoundary == b2.uBoundary && b1.pBoundary == b2.pBoundary;
}
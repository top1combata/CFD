#pragma once

#include <Eigen/Core>
#include <vector>
#include <array>

using Index = int;

template<class T>
using List = std::vector<T>;

template<class T, Index n>
using Array = std::array<T,n>;


using Scalar = double;

using Vector = Eigen::Matrix<Scalar, 3, 1>;

using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using ScalarField = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

using VectorField = Eigen::Matrix<Scalar, Eigen::Dynamic, 3>;
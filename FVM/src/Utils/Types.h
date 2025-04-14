#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include <array>
#include <unordered_map>


using Index = int;

template<class T>
using List = std::vector<T>;

template<class T, Index n>
using Array = std::array<T,n>;

template<class K, class V, class Hash = std::hash<K>, class Equal = std::equal_to<K>>
using HashMap = std::unordered_map<K, V, Hash, Equal>;



using Scalar = double;

using Vector = Eigen::Matrix<Scalar, 3, 1>;

using Tensor = Eigen::Matrix<Scalar, 3, 3>;

using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;

template<class T>
using Field = Eigen::Matrix<T, Eigen::Dynamic, 1>;

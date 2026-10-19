#pragma once
#include "euler/types.hpp"
namespace euler {
// Allocation-free small GEMM, with contiguous output rows handled by Eigen.
// Destinations must be sized and must not alias either input.
inline void multiply(ConstBlockRef left, ConstBlockRef right, BlockRef result) {
    for (int j = 0; j < right.cols(); ++j) {
        result.col(j).setZero();
        for (int k = 0; k < left.cols(); ++k)
            result.col(j).noalias() += right(k, j) * left.col(k);
    }
}
inline void multiply_transpose_right(ConstBlockRef left, ConstBlockRef right, BlockRef result) {
    for (int j = 0; j < right.rows(); ++j) {
        result.col(j).setZero();
        for (int k = 0; k < left.cols(); ++k)
            result.col(j).noalias() += right(j, k) * left.col(k);
    }
}
} // namespace euler

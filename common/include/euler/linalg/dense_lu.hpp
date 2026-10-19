#pragma once
#include "euler/types.hpp"
#include <limits>
namespace euler {
// Eigen owns the pivoted factorization. The wrapper only manages preallocated
// RHS scratch and alias-safe vector/right solves; no custom LU algorithm.
class DenseLU {
  public:
    explicit DenseLU(int width) : lu_(width), scratch_(width) {}
    void compute(ConstBlockRef matrix) {
        require(matrix.allFinite(), "Nonfinite block LU input");
        lu_.compute(matrix);
        require(lu_.matrixLU().allFinite() &&
                    (lu_.matrixLU().diagonal().array().abs() > std::numeric_limits<double>::min()).all(),
                "Singular or nonfinite LU factors");
    }
    template <class Input, class Output> void solve(const Input &b, Output &&x) const {
        scratch_.noalias() = lu_.permutationP() * b;
        lu_.matrixLU().template triangularView<Eigen::UnitLower>().solveInPlace(scratch_);
        lu_.matrixLU().template triangularView<Eigen::Upper>().solveInPlace(scratch_);
        x = scratch_;
    }
    // In place right solve: X <- X A^-1, one row at a time using U^T,L^T,P.
    void right_solve(BlockRef x) const {
        for (int row = 0; row < x.rows(); ++row) {
            scratch_ = x.row(row).transpose();
            lu_.matrixLU().transpose().template triangularView<Eigen::Lower>().solveInPlace(scratch_);
            lu_.matrixLU().transpose().template triangularView<Eigen::UnitUpper>().solveInPlace(scratch_);
            x.row(row).transpose().noalias() = lu_.permutationP().transpose() * scratch_;
        }
    }

  private:
    Eigen::PartialPivLU<Mat> lu_;
    mutable Vec scratch_;
};
// Small dense block update without Eigen GEMM's dynamically allocated packing.
inline void subtract_product(BlockRef destination, ConstBlockRef left, ConstBlockRef right) {
    for (int j = 0; j < right.cols(); ++j)
        for (int k = 0; k < left.cols(); ++k)
            destination.col(j).noalias() -= right(k, j) * left.col(k);
}
} // namespace euler

#pragma once
#include "block_graph.hpp"
#include "dense_lu.hpp"
#include "gmres.hpp"
namespace euler {
// Small-problem reference. Expand block connectivity into a preallocated dense
// column-major tensor [row_node,row_variable,row_cell,col_node,col_variable,col_cell].
// Full scalar partial pivoting is retained, including pivots across cell blocks.
// O(ndof^2) storage / O(ndof^3) work is deliberate for this reference path.
class DirectSolver {
  public:
    DirectSolver(const BlockGraph &graph, int width)
        : graph_(graph), width_(width), matrix_(graph.rows * width, graph.rows * width),
          block_(width, width), lu_(graph.rows * width), residual_(graph.rows * width) {}
    template <class BlockProvider> void refresh(BlockProvider &&provider) {
        matrix_.setZero();
        for (int row = 0; row < graph_.rows; ++row)
            for (int p = graph_.row_ptr[row]; p < graph_.row_ptr[row + 1]; ++p) {
                provider(row, graph_.col_idx[p], block_);
                matrix_.block(row * width_, graph_.col_idx[p] * width_, width_, width_) = block_;
            }
        lu_.compute(matrix_);
    }
    LinearResult solve(const Vec &rhs, Vec &x, const Vec &inverse_mass) {
        lu_.solve(rhs, x);
        residual_.noalias() = matrix_ * x;
        residual_ -= rhs;
        double denominator = (rhs.array().square() * inverse_mass.array()).sum();
        double numerator = (residual_.array().square() * inverse_mass.array()).sum();
        double relative =
            denominator > 0 ? std::sqrt(numerator / denominator) : std::sqrt(numerator);
        require(std::isfinite(relative), "Nonfinite direct residual");
        return {0, relative};
    }

  private:
    BlockGraph graph_;
    int width_;
    Mat matrix_, block_;
    DenseLU lu_;
    Vec residual_;
};
} // namespace euler

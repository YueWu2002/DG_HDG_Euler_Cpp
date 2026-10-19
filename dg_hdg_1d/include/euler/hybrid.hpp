#pragma once
#include "volume_jacobian.hpp"
#include "linalg/dense_lu.hpp"
namespace euler {
struct CellCondensation {
    Mat A, B, C, X; // X=A^-1 B by solves; never store A^-1.
    Vec y;
    DenseLU lu;
    explicit CellCondensation(int width)
        : A(width, width), B(width, NFACES * NVAR), C(NFACES * NVAR, width),
          X(width, NFACES * NVAR), y(width), lu(width) {}
};
class HybridOperator {
  public:
    explicit HybridOperator(const Problem &problem);
    void residual(const Vec &u, const Vec &lambda, Vec &volume, Vec &skeleton) const;
    void refresh(const Vec &u, const Vec &lambda);
    void apply(const Vec &x, Vec &result) const;
    void block(int row, int col, BlockRef result) const;
    void recover(const Vec &trace_direction, Vec &volume_direction) const;
    const Vec &rhs() const { return rhs_; }

  private:
    const Problem &p_;
    VolumeJacobian volume_;
    std::vector<CellCondensation> cells_;
    Vec volume_residual_, skeleton_residual_, rhs_;
    mutable Vec local_, faces_, face_result_;
    int face_index(int cell, int side) const;
    void gather(const Vec &x, int cell) const;
};
} // namespace euler

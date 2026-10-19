#pragma once
#include "volume_jacobian.hpp"
#include "face_jacobian.hpp"
#include "euler/linalg/dense_lu.hpp"
namespace euler {
class HybridOperator {
  public:
    explicit HybridOperator(const Problem &problem);
    void residual(const Vec &state, const Vec &lambda, Vec &volume, Vec &skeleton) const;
    void refresh(const Vec &state, const Vec &lambda);
    void apply(const Vec &direction, Vec &result) const;
    void block(int row, int col, BlockRef result) const;
    void recover(const Vec &trace_direction, Vec &volume_direction) const;
    const Vec &rhs() const { return rhs_; }

  private:
    const Problem &p_;
    VolumeJacobian volume_;
    FaceJacobian faces_jacobian_;
    std::vector<DenseLU> lu_;
    // Contiguous cell-blocked pools: [row_dof,col_dof,cell].
    Mat x_, schur_, y_, a_, b_;
    Vec volume_residual_, skeleton_residual_, rhs_;
    mutable Vec local_, faces_, face_result_;
    void gather(const Vec &input, int cell) const;
    void scatter(Vec &output, int cell, const Vec &local) const;
    void assemble_cell(int cell);
};
} // namespace euler

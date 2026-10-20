#pragma once
#include "volume_jacobian.hpp"
#include "face_jacobian.hpp"
#include "euler/linalg/dense_lu.hpp"
namespace euler {
class HybridOperator {
  public:
    explicit HybridOperator(const Problem &problem);
    void residual(const Vec &state,const Vec &lambda,Vec &volume,Vec &skeleton) const;
    void refresh(const Vec &state,const Vec &lambda);
    void refresh_projected(const Vec &state,const Vec &dg_residual);
    void apply(const Vec &direction,Vec &result) const;
    void block(int row,int col,BlockRef result) const;
    void recover(const Vec &trace_direction,Vec &volume_direction) const;
    void recovered_residual_from_schur(const Vec &equation_error,Vec &result) const;
    const Vec &rhs() const { return rhs_; }
  private:
    const Problem &p_;
    VolumeJacobian volume_;
    FaceJacobian faces_jacobian_;
    DenseLU lu_;
    // Only X, local Schur and y survive a cell. A/LU/RHS reuse one workspace.
    Vec x_, schur_;
    Mat y_, a_, local_rhs_, solved_, face_action_;
    Vec volume_residual_, skeleton_residual_, rhs_;
    mutable Vec local_, faces_, face_result_;
    mutable Mat face_flux_;
    void gather(const Vec &input,int cell) const;
    void scatter(Vec &output,int cell,const Vec &local) const;
    void assemble_cell(int cell);
    void condense(const Vec &state);
};
}

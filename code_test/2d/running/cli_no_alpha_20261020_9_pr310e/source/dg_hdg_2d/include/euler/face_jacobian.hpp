#pragma once
#include "discretization.hpp"
namespace euler {
class FaceJacobian {
  public:
    explicit FaceJacobian(const Problem &problem);
    void refresh(const Vec &state);
    void hll_action(int face, ConstBlockRef left, ConstBlockRef right, BlockRef flux) const;
    FluxJacobian hll(int face, int q, bool left) const;
    FluxJacobian hdg(int face, int q, bool left) const;
    double tau(int face, int q, bool left) const;
    void cell_c(int cell, Eigen::Ref<const Vec> direction, Eigen::Ref<Vec> result) const;
    void cell_c_many(int cell, ConstBlockRef directions, BlockRef result) const;
    const Vec &trace_diagonal() const { return diagonal_; }
  private:
    const Problem &p_;
    // Tangent node contiguous; cache u,v,H instead of two full 4x4 derivatives.
    Mat primitives_left_, primitives_right_, tau_left_, tau_right_, weights_;
    Mat blend_left_, blend_right_, penalty_;
    Vec diagonal_;
    mutable Mat left_, right_, flux_left_, flux_right_, samples_, coefficients_;
    mutable Vec pressure_direction_;
    void normal_action(const Mat &primitive,int face,ConstBlockRef direction,BlockRef result,int axis) const;
    FluxJacobian physical(int face,int q,bool left) const;
};
}

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
    const Vec &trace_diagonal() const { return diagonal_; }

  private:
    const Problem &p_;
    // [tangent_node,flux_variable,state_variable,face].
    Mat physical_left_, physical_right_, hll_left_, hll_right_, tau_left_, tau_right_, weights_;
    Vec diagonal_;
    mutable Mat left_, right_;
};
} // namespace euler

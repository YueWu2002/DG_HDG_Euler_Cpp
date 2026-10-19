#pragma once
#include "discretization.hpp"
namespace euler {
// Partial assembly: store only F'(U_q), never a cell-sized dense Jacobian.
// refresh() snapshots one Newton state; apply() and block() reuse that snapshot.
// All storage, including single-cell scratch, is acquired in the constructor.
class VolumeJacobian {
  public:
    explicit VolumeJacobian(const Problem &problem);
    void refresh(const Vec &state);
    void apply(const Vec &direction, Vec &result) const;
    void block(int cell, BlockRef result) const;

  private:
    const Problem &p_;
    // Column-major tensor [q,flux_variable,state_variable,cell].
    // Each coefficient is contiguous across quadrature points for SIMD.
    Mat coefficients_, derivative_transpose_;
    Mat point_states_;
    mutable Mat point_direction_, flux_direction_;
    mutable Vec point_weights_;
    double inverse_step_ = 0;
    bool ready_ = false;
    void store_coefficients(int cell, ConstBlockRef point_states);
    void add_cell_action(int cell, ConstBlockRef point_direction, BlockRef result) const;
    template <int Nodes, int Points>
    void add_small_cell_action(int cell, ConstBlockRef point_direction, BlockRef result) const;
    template <int Nodes, int Points>
    void add_overintegrated_cell_action(int cell, ConstBlockRef nodal, BlockRef result) const;
};
} // namespace euler

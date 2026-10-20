#pragma once
#include "discretization.hpp"
namespace euler {
class VolumeJacobian {
  public:
    explicit VolumeJacobian(const Problem &problem);
    void refresh(const Vec &state);
    void apply(const Vec &direction, Vec &result) const;
    void block(int cell, BlockRef result) const;

  private:
    const Problem &p_;
    TensorProduct tensor_;
    // [qx,qy,{u,v,H},cell]: 3 cached fields instead of 32 Jacobian entries.
    Mat primitives_, derivative_transpose_, weighted_values_transpose_;
    mutable Mat x_, y_; // One cell of explicit derivatives, only for block queries.
    mutable Mat points_, point_direction_, flux_x_, flux_y_, line_x_, line_y_;
    mutable Vec weights_x_, weights_y_, pressure_direction_;
    double inverse_step_ = 0;
    bool ready_ = false;
    void store(int cell, ConstBlockRef states);
    void cache_cell_jacobians(int cell) const;
    void point_action(int cell, ConstBlockRef direction) const;
    void collocated_block(int cell, BlockRef result) const;
    void overintegrated_block(int cell, BlockRef result) const;
};
} // namespace euler

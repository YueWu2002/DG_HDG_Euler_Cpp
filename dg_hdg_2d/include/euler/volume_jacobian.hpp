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
    // [qx,qy,flux_variable,state_variable,cell], qx contiguous.
    Mat x_, y_, derivative_transpose_, weighted_values_transpose_;
    mutable Mat points_, point_direction_, flux_x_, flux_y_, line_x_, line_y_;
    mutable Vec weights_x_, weights_y_;
    double inverse_step_ = 0;
    bool ready_ = false;
    void store(int cell, ConstBlockRef states);
    void point_action(int cell, ConstBlockRef direction) const;
    void collocated_block(int cell, BlockRef result) const;
    void overintegrated_block(int cell, BlockRef result) const;
};
} // namespace euler

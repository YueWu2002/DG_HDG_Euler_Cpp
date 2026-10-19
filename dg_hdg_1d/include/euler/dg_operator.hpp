#pragma once
#include "volume_jacobian.hpp"
namespace euler {
// Snapshot of the DG Newton Jacobian, shared by Jv and block assembly.
// Construct before stepping; refresh after setting dt, frozen speeds, and state.
class DGOperator {
  public:
    explicit DGOperator(const Problem &problem);
    void refresh(const Vec &state);
    void apply(const Vec &direction, Vec &result) const;
    void block(int row, int col, BlockRef result) const;

  private:
    const Problem &p_;
    VolumeJacobian volume_;
    // Column-major tensor [flux_variable,state_variable,face].
    Mat face_left_, face_right_;
    bool ready_ = false;
};
} // namespace euler

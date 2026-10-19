#pragma once
#include "volume_jacobian.hpp"
#include "face_jacobian.hpp"
namespace euler {
class DGOperator {
  public:
    explicit DGOperator(const Problem &problem);
    void refresh(const Vec &state);
    void apply(const Vec &direction, Vec &result) const;
    void block(int row, int col, BlockRef result) const;

  private:
    const Problem &p_;
    VolumeJacobian volume_;
    FaceJacobian faces_;
    mutable Mat left_, right_, flux_;
};
} // namespace euler

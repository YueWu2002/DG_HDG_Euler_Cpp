#pragma once
#include "euler/nodal_basis.hpp"
namespace euler {
// Sum factorization on [x,y] sheets. No n^2 by Q^2 interpolation matrix.
class TensorProduct {
  public:
    TensorProduct(const NodalBasis &basis, double hx, double hy);
    void interpolate(ConstBlockRef nodal, BlockRef points) const;
    void subtract_weak(ConstBlockRef flux_x, ConstBlockRef flux_y, BlockRef residual) const;

  private:
    template<int N> void subtract_collocated(ConstBlockRef fx, ConstBlockRef fy, BlockRef residual) const;
    const NodalBasis &b_;
    double scale_x_, scale_y_;
    Mat derivative_transpose_, weighted_values_, weighted_values_transpose_;
    mutable Mat interpolation_, sheet_, scratch_, weak_;
};
} // namespace euler

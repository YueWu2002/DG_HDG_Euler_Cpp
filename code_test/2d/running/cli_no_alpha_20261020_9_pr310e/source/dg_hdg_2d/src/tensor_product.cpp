#include "euler/tensor_product.hpp"
#include "euler/tensor_kernels.hpp"
namespace euler {
TensorProduct::TensorProduct(const NodalBasis &b, double hx, double hy)
    : b_(b), scale_x_(.5 * hy), scale_y_(.5 * hx),
      derivative_transpose_(b.weighted_derivatives.transpose()),
      weighted_values_(b.qw.asDiagonal() * b.values),
      weighted_values_transpose_(weighted_values_.transpose()),
      interpolation_(b.points, b.nodes), sheet_(b.points, b.points),
      scratch_(b.nodes, b.points), weak_(b.nodes, b.nodes) {}
void TensorProduct::interpolate(ConstBlockRef nodal, BlockRef points) const {
    if (b_.points == b_.nodes) {
        points = nodal;
        return;
    }
    for (int v = 0; v < NVAR; ++v) {
        Eigen::Map<const Mat> input(nodal.col(v).data(), b_.nodes, b_.nodes);
        multiply(b_.values, input, interpolation_);
        multiply_transpose_right(interpolation_, b_.values, sheet_);
        points.col(v) = Eigen::Map<const Vec>(sheet_.data(), b_.points * b_.points);
    }
}
template<int N>
void TensorProduct::subtract_collocated(ConstBlockRef fx, ConstBlockRef fy, BlockRef residual) const {
    using Sheet=Eigen::Matrix<double,N,N>;
    const Eigen::Map<const Sheet> d(b_.weighted_derivatives.data());
    for(int v=0;v<NVAR;++v) {
        const Eigen::Map<const Sheet> x(fx.col(v).data()), y(fy.col(v).data());
        Eigen::Map<Sheet> output(residual.col(v).data());
        Sheet wx,wy;
        wx.noalias()=d.transpose()*x;
        wy.noalias()=y*d;
        for(int j=0;j<N;++j) {
            output.col(j).noalias()-=(scale_x_*b_.w[j])*wx.col(j);
            output.col(j).array()-=scale_y_*b_.w.array()*wy.col(j).array();
        }
    }
}
void TensorProduct::subtract_weak(ConstBlockRef fx, ConstBlockRef fy, BlockRef residual) const {
    if (b_.nodes == 1)
        return;
    if(b_.points==b_.nodes) {
        switch(b_.nodes) {
        case 2: subtract_collocated<2>(fx,fy,residual); return;
        case 3: subtract_collocated<3>(fx,fy,residual); return;
        case 4: subtract_collocated<4>(fx,fy,residual); return;
        case 5: subtract_collocated<5>(fx,fy,residual); return;
        case 6: subtract_collocated<6>(fx,fy,residual); return;
        case 7: subtract_collocated<7>(fx,fy,residual); return;
        case 8: subtract_collocated<8>(fx,fy,residual); return;
        case 9: subtract_collocated<9>(fx,fy,residual); return;
        }
    }
    for (int v = 0; v < NVAR; ++v) {
        Eigen::Map<const Mat> x(fx.col(v).data(), b_.points, b_.points);
        Eigen::Map<const Mat> y(fy.col(v).data(), b_.points, b_.points);
        Eigen::Map<Mat> output(residual.col(v).data(), b_.nodes, b_.nodes);
        multiply(derivative_transpose_, x, scratch_);
        if (b_.points == b_.nodes) {
            for (int j = 0; j < b_.nodes; ++j)
                output.col(j).noalias() -= (scale_x_ * b_.w[j]) * scratch_.col(j);
            multiply(y, b_.weighted_derivatives, weak_);
            for (int j = 0; j < b_.nodes; ++j)
                output.col(j).array() -= scale_y_ * b_.w.array() * weak_.col(j).array();
        } else {
            multiply(scratch_, weighted_values_, weak_);
            output.noalias() -= scale_x_ * weak_;
            multiply(weighted_values_transpose_, y, scratch_);
            multiply(scratch_, b_.weighted_derivatives, weak_);
            output.noalias() -= scale_y_ * weak_;
        }
    }
}
} // namespace euler

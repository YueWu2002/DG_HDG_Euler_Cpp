#include "euler/volume_jacobian.hpp"
namespace euler {
VolumeJacobian::VolumeJacobian(const Problem &p)
    : p_(p), tensor_(p.basis, p.hx, p.hy),
      primitives_(p.basis.points*p.basis.points, 3*p.layout.cells),
      derivative_transpose_(p.basis.weighted_derivatives.transpose()),
      weighted_values_transpose_(p.basis.values.transpose()*p.basis.qw.asDiagonal()),
      x_(primitives_.rows(),NVAR*NVAR), y_(x_.rows(),x_.cols()),
      points_(x_.rows(), NVAR), point_direction_(x_.rows(), NVAR),
      flux_x_(x_.rows(), NVAR), flux_y_(x_.rows(), NVAR),
      line_x_(p.basis.nodes, p.basis.nodes), line_y_(p.basis.nodes, p.basis.nodes),
      weights_x_(p.basis.points), weights_y_(p.basis.points), pressure_direction_(x_.rows()) {}
void VolumeJacobian::store(int cell, ConstBlockRef states) {
    for(int q=0;q<states.rows();++q) {
        const auto primitive=p_.example.gas.primitive(states.row(q).transpose());
        primitives_(q,3*cell)=primitive.velocity[0];
        primitives_(q,3*cell+1)=primitive.velocity[1];
        primitives_(q,3*cell+2)=primitive.enthalpy;
    }
}
void VolumeJacobian::cache_cell_jacobians(int cell) const {
    for(int q=0;q<primitives_.rows();++q) {
        const auto x=p_.example.gas.cached_jacobian(primitives_(q,3*cell),primitives_(q,3*cell+1),primitives_(q,3*cell+2),0);
        const auto y=p_.example.gas.cached_jacobian(primitives_(q,3*cell),primitives_(q,3*cell+1),primitives_(q,3*cell+2),1);
        for(int w=0;w<NVAR;++w) for(int v=0;v<NVAR;++v) {
            x_(q,w*NVAR+v)=x(v,w); y_(q,w*NVAR+v)=y(v,w);
        }
    }
}
void VolumeJacobian::refresh(const Vec &state) {
    require(p_.step_size > 0, "Set the time step before refreshing the Jacobian");
    inverse_step_ = 1 / p_.step_size;
    if (p_.basis.nodes > 1)
        for (int cell = 0; cell < p_.layout.cells; ++cell) {
            if (p_.basis.points == p_.basis.nodes) store(cell, p_.layout.cell(state, cell));
            else {
                tensor_.interpolate(p_.layout.cell(state, cell), points_);
                store(cell, points_);
            }
        }
    ready_ = true;
}
void VolumeJacobian::point_action(int cell, ConstBlockRef direction) const {
    const auto u=primitives_.col(3*cell).array(), v=primitives_.col(3*cell+1).array();
    const auto H=primitives_.col(3*cell+2).array();
    const auto dr=direction.col(Density).array(), dx=direction.col(MomentumX).array();
    const auto dy=direction.col(MomentumY).array(), de=direction.col(Energy).array();
    pressure_direction_.array()=(p_.example.gas.gamma()-1.)*(de-u*dx-v*dy+.5*(u.square()+v.square())*dr);
    const auto dp=pressure_direction_.array();
    flux_x_.col(Density)=direction.col(MomentumX);
    flux_y_.col(Density)=direction.col(MomentumY);
    flux_x_.col(MomentumX).array()=2*u*dx-u.square()*dr+dp;
    flux_x_.col(MomentumY).array()=u*dy+v*dx-u*v*dr;
    flux_y_.col(MomentumX)=flux_x_.col(MomentumY);
    flux_y_.col(MomentumY).array()=2*v*dy-v.square()*dr+dp;
    flux_x_.col(Energy).array()=H*(dx-u*dr)+u*(de+dp);
    flux_y_.col(Energy).array()=H*(dy-v*dr)+v*(de+dp);
}
void VolumeJacobian::apply(const Vec &direction, Vec &result) const {
    require(ready_, "Refresh the volume Jacobian before applying it");
    result.array() = p_.mass.array() * direction.array() * inverse_step_;
    if (p_.basis.nodes == 1) return;
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        if (p_.basis.points == p_.basis.nodes) point_action(cell, p_.layout.cell(direction, cell));
        else {
            tensor_.interpolate(p_.layout.cell(direction, cell), point_direction_);
            point_action(cell, point_direction_);
        }
        tensor_.subtract_weak(flux_x_, flux_y_, p_.layout.cell(result, cell));
    }
}
void VolumeJacobian::collocated_block(int /*cell*/, BlockRef result) const {
    const int n = p_.basis.nodes, plane = p_.layout.plane();
    for (int w = 0; w < NVAR; ++w)
        for (int jy = 0; jy < n; ++jy)
            for (int jx = 0; jx < n; ++jx) {
                int node = jx + n * jy;
                auto column = result.col(p_.layout.local(w, node));
                for (int v = 0; v < NVAR; ++v) {
                    int coefficient = w * NVAR + v;
                    column.segment(v * plane + n * jy, n).noalias() -=
                        (.5 * p_.hy * p_.basis.w[jy] * x_(node, coefficient)) *
                        derivative_transpose_.col(jx);
                    double scale = -.5 * p_.hx * p_.basis.w[jx] * y_(node, coefficient);
                    for (int iy = 0; iy < n; ++iy)
                        column[v * plane + jx + n * iy] += scale * derivative_transpose_(iy, jy);
                }
            }
}
void VolumeJacobian::overintegrated_block(int /*cell*/, BlockRef result) const {
    const auto &b = p_.basis;
    const int n = b.nodes, Q = b.points, plane = p_.layout.plane();
    // Contract qx first into n-by-n line matrices; then assemble qy factors.
    // This avoids a Q^2-by-n^4 quadrature expansion.
    for (int qb = 0; qb < Q; ++qb)
        for (int w = 0; w < NVAR; ++w)
            for (int v = 0; v < NVAR; ++v) {
                int coefficient = w * NVAR + v;
                for (int k = 0; k < n; ++k) {
                    weights_x_.array() = x_.col(coefficient).segment(qb * Q, Q).array() * b.values.col(k).array();
                    weights_y_.array() = y_.col(coefficient).segment(qb * Q, Q).array() * b.values.col(k).array();
                    line_x_.col(k).noalias() = -derivative_transpose_ * weights_x_;
                    line_y_.col(k).noalias() = -weighted_values_transpose_ * weights_y_;
                }
                for (int l = 0; l < n; ++l)
                    for (int k = 0; k < n; ++k) {
                        auto column = result.col(p_.layout.local(w, k + n * l));
                        for (int j = 0; j < n; ++j) {
                            double cx = .5 * p_.hy * b.qw[qb] * b.values(qb, j) * b.values(qb, l);
                            double cy = .5 * p_.hx * b.weighted_derivatives(qb, j) * b.values(qb, l);
                            column.segment(v * plane + n * j, n).noalias() +=
                                cx * line_x_.col(k) + cy * line_y_.col(k);
                        }
                    }
            }
}
void VolumeJacobian::block(int cell, BlockRef result) const {
    require(ready_, "Refresh the volume Jacobian before assembling a block");
    result.setZero();
    if (p_.basis.nodes > 1) {
        cache_cell_jacobians(cell);
        if (p_.basis.points == p_.basis.nodes) collocated_block(cell, result);
        else overintegrated_block(cell, result);
    }
    result.diagonal().array() += p_.mass.segment(cell * p_.layout.block_size(), p_.layout.block_size()).array() * inverse_step_;
}
} // namespace euler

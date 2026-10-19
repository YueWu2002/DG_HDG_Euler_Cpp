#include "euler/volume_jacobian.hpp"
namespace euler {
VolumeJacobian::VolumeJacobian(const Problem &p)
    : p_(p), tensor_(p.basis, p.hx, p.hy),
      x_(p.basis.points * p.basis.points, NVAR * NVAR * p.layout.cells), y_(x_.rows(), x_.cols()),
      derivative_transpose_(p.basis.weighted_derivatives.transpose()),
      weighted_values_transpose_(p.basis.values.transpose() * p.basis.qw.asDiagonal()),
      points_(x_.rows(), NVAR), point_direction_(x_.rows(), NVAR),
      flux_x_(x_.rows(), NVAR), flux_y_(x_.rows(), NVAR),
      line_x_(p.basis.nodes, p.basis.nodes), line_y_(p.basis.nodes, p.basis.nodes),
      weights_x_(p.basis.points), weights_y_(p.basis.points) {}
void VolumeJacobian::store(int cell, ConstBlockRef states) {
    for (int q = 0; q < states.rows(); ++q) {
        State state = states.row(q).transpose();
        auto primitive = p_.example.gas.primitive(state);
        FluxJacobian x = p_.example.gas.jacobian(state, primitive, 0);
        FluxJacobian y = p_.example.gas.jacobian(state, primitive, 1);
        for (int w = 0; w < NVAR; ++w)
            for (int v = 0; v < NVAR; ++v) {
                int column = (cell * NVAR + w) * NVAR + v;
                x_(q, column) = x(v, w);
                y_(q, column) = y(v, w);
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
    for (int v = 0; v < NVAR; ++v) {
        flux_x_.col(v).setZero();
        flux_y_.col(v).setZero();
        for (int w = 0; w < NVAR; ++w) {
            int column = (cell * NVAR + w) * NVAR + v;
            flux_x_.col(v).array() += x_.col(column).array() * direction.col(w).array();
            flux_y_.col(v).array() += y_.col(column).array() * direction.col(w).array();
        }
    }
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
void VolumeJacobian::collocated_block(int cell, BlockRef result) const {
    const int n = p_.basis.nodes, plane = p_.layout.plane();
    for (int w = 0; w < NVAR; ++w)
        for (int jy = 0; jy < n; ++jy)
            for (int jx = 0; jx < n; ++jx) {
                int node = jx + n * jy;
                auto column = result.col(p_.layout.local(w, node));
                for (int v = 0; v < NVAR; ++v) {
                    int coefficient = (cell * NVAR + w) * NVAR + v;
                    column.segment(v * plane + n * jy, n).noalias() -=
                        (.5 * p_.hy * p_.basis.w[jy] * x_(node, coefficient)) *
                        derivative_transpose_.col(jx);
                    double scale = -.5 * p_.hx * p_.basis.w[jx] * y_(node, coefficient);
                    for (int iy = 0; iy < n; ++iy)
                        column[v * plane + jx + n * iy] += scale * derivative_transpose_(iy, jy);
                }
            }
}
void VolumeJacobian::overintegrated_block(int cell, BlockRef result) const {
    const auto &b = p_.basis;
    const int n = b.nodes, Q = b.points, plane = p_.layout.plane();
    // Contract qx first into n-by-n line matrices; then assemble qy factors.
    // This avoids a Q^2-by-n^4 quadrature expansion.
    for (int qb = 0; qb < Q; ++qb)
        for (int w = 0; w < NVAR; ++w)
            for (int v = 0; v < NVAR; ++v) {
                int coefficient = (cell * NVAR + w) * NVAR + v;
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
        if (p_.basis.points == p_.basis.nodes) collocated_block(cell, result);
        else overintegrated_block(cell, result);
    }
    result.diagonal().array() += p_.mass.segment(cell * p_.layout.block_size(), p_.layout.block_size()).array() * inverse_step_;
}
} // namespace euler

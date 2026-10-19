#include "euler/volume_jacobian.hpp"
namespace euler {
VolumeJacobian::VolumeJacobian(const Problem &p)
    : p_(p), coefficients_(p.basis.points, NVAR * NVAR * p.layout.cells),
      derivative_transpose_(p.basis.weighted_derivatives.transpose()),
      point_states_(p.basis.points, NVAR), point_direction_(p.basis.points, NVAR),
      flux_direction_(p.basis.points, NVAR), point_weights_(p.basis.points) {}

void VolumeJacobian::store_coefficients(int cell, ConstBlockRef states) {
    for (int q = 0; q < p_.basis.points; ++q) {
        FluxJacobian jacobian = p_.example.gas.jacobian(states.row(q).transpose());
        for (int w = 0; w < NVAR; ++w)
            for (int v = 0; v < NVAR; ++v)
                coefficients_(q, (cell * NVAR + w) * NVAR + v) = jacobian(v, w);
    }
}
void VolumeJacobian::refresh(const Vec &state) {
    require(p_.step_size > 0, "Set the time step before refreshing the volume Jacobian");
    inverse_step_ = 1 / p_.step_size;
    // A constant test basis has zero derivative, even with overintegration.
    if (p_.layout.nodes == 1) {
        ready_ = true;
        return;
    }
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        auto nodal = p_.layout.cell(state, cell);
        if (p_.basis.points == p_.layout.nodes) {
            // Collocation: V=I. No interpolation or copy of nodal states.
            store_coefficients(cell, nodal);
        } else {
            for (int v = 0; v < NVAR; ++v)
                point_states_.col(v).noalias() = p_.basis.values * nodal.col(v);
            store_coefficients(cell, point_states_);
        }
    }
    ready_ = true;
}
void VolumeJacobian::add_cell_action(int cell, ConstBlockRef direction, BlockRef result) const {
    for (int v = 0; v < NVAR; ++v) {
        flux_direction_.col(v).setZero();
        for (int w = 0; w < NVAR; ++w)
            flux_direction_.col(v).array() +=
                coefficients_.col((cell * NVAR + w) * NVAR + v).array() *
                direction.col(w).array();
        result.col(v).noalias() -= derivative_transpose_ * flux_direction_.col(v);
    }
}
template <int Nodes, int Points>
void VolumeJacobian::add_small_cell_action(int cell, ConstBlockRef direction,
                                          BlockRef result) const {
    // Small fixed-size Eigen kernels avoid dynamic-loop overhead and expose
    // the contiguous quadrature axis to SIMD. All maps here are contiguous.
    using Cell = Eigen::Matrix<double, Nodes, NVAR, Eigen::ColMajor>;
    using PointStates = Eigen::Matrix<double, Points, NVAR, Eigen::ColMajor>;
    using Derivative = Eigen::Matrix<double, Nodes, Points, Eigen::ColMajor>;
    Eigen::Map<const PointStates> input(direction.data());
    Eigen::Map<Cell> output(result.data());
    Eigen::Map<const Derivative> derivative(derivative_transpose_.data());
    PointStates flux = PointStates::Zero();
    for (int v = 0; v < NVAR; ++v)
        for (int w = 0; w < NVAR; ++w) {
            Eigen::Map<const Eigen::Array<double, Points, 1>> coefficient(
                coefficients_.col((cell * NVAR + w) * NVAR + v).data());
            flux.col(v).array() += coefficient * input.col(w).array();
        }
    output.noalias() -= derivative * flux;
}
template <int Nodes, int Points>
void VolumeJacobian::add_overintegrated_cell_action(int cell, ConstBlockRef nodal,
                                                   BlockRef result) const {
    using Interpolation = Eigen::Matrix<double, Points, Nodes, Eigen::ColMajor>;
    using Cell = Eigen::Matrix<double, Nodes, NVAR, Eigen::ColMajor>;
    Eigen::Map<const Interpolation> interpolation(p_.basis.values.data());
    Eigen::Map<const Cell> input(nodal.data());
    Eigen::Matrix<double, Points, NVAR> direction = interpolation * input;
    add_small_cell_action<Nodes, Points>(cell, direction, result);
}
void VolumeJacobian::apply(const Vec &direction, Vec &result) const {
    require(ready_, "Refresh the volume Jacobian before applying it");
    result.array() = p_.mass.array() * direction.array() * inverse_step_;
    if (p_.layout.nodes == 1)
        return;
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        auto nodal = p_.layout.cell(direction, cell);
        auto output = p_.layout.cell(result, cell);
        if (p_.basis.points == p_.layout.nodes) {
            // Common P1/P3 kernels; all other degrees retain the general path.
            if (p_.layout.nodes == 2)
                add_small_cell_action<2, 2>(cell, nodal, output);
            else if (p_.layout.nodes == 4)
                add_small_cell_action<4, 4>(cell, nodal, output);
            else
                add_cell_action(cell, nodal, output);
        } else if (p_.layout.nodes == 2 && p_.basis.points == 4) {
            add_overintegrated_cell_action<2, 4>(cell, nodal, output);
        } else if (p_.layout.nodes == 4 && p_.basis.points == 6) {
            add_overintegrated_cell_action<4, 6>(cell, nodal, output);
        } else {
            for (int v = 0; v < NVAR; ++v)
                point_direction_.col(v).noalias() = p_.basis.values * nodal.col(v);
            add_cell_action(cell, point_direction_, output);
        }
    }
}
void VolumeJacobian::block(int cell, BlockRef result) const {
    require(ready_, "Refresh the volume Jacobian before assembling a block");
    const int nodes = p_.layout.nodes;
    if (nodes == 1) {
        result.setZero();
        result.diagonal().array() =
            p_.mass.segment(cell * NVAR, NVAR).array() * inverse_step_;
        return;
    }
    // Column first, then contiguous row nodes. Every entry is overwritten once.
    for (int w = 0; w < NVAR; ++w)
        for (int j = 0; j < nodes; ++j)
            for (int v = 0; v < NVAR; ++v) {
                auto output = result.col(p_.layout.local(w, j)).segment(v * nodes, nodes);
                auto coefficient = coefficients_.col((cell * NVAR + w) * NVAR + v);
                if (p_.basis.points == nodes) {
                    // A_(v,i),(w,j) = -D^w_(j,i) J_j(v,w).
                    output.noalias() = -coefficient[j] * derivative_transpose_.col(j);
                } else {
                    point_weights_.array() = coefficient.array() * p_.basis.values.col(j).array();
                    output.noalias() = -derivative_transpose_ * point_weights_;
                }
            }
    result.diagonal().array() +=
        p_.mass.segment(cell * p_.layout.block_size(), p_.layout.block_size()).array() *
        inverse_step_;
}
} // namespace euler

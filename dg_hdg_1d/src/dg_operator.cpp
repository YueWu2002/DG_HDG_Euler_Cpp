#include "euler/dg_operator.hpp"
namespace euler {
DGOperator::DGOperator(const Problem &p)
    : p_(p), volume_(p), face_left_(NVAR * NVAR, p.layout.cells),
      face_right_(NVAR * NVAR, p.layout.cells) {}

void DGOperator::refresh(const Vec &state) {
    volume_.refresh(state);
    // Each face is evaluated once, independently of the block graph/ordering.
    for (int face = 0; face < p_.layout.cells; ++face) {
        int next = (face + 1) % p_.layout.cells;
        auto speed = p_.speeds[face];
        Eigen::Map<FluxJacobian> left(face_left_.col(face).data());
        Eigen::Map<FluxJacobian> right(face_right_.col(face).data());
        left = speed.b / (speed.a + speed.b) *
               (p_.example.gas.jacobian(p_.sample(state, face, p_.basis.right)) +
                speed.a * FluxJacobian::Identity());
        right = speed.a / (speed.a + speed.b) *
                (p_.example.gas.jacobian(p_.sample(state, next, p_.basis.left)) -
                 speed.b * FluxJacobian::Identity());
    }
    ready_ = true;
}
void DGOperator::apply(const Vec &direction, Vec &result) const {
    require(ready_, "Refresh the DG Jacobian before applying it");
    volume_.apply(direction, result);
    if (p_.layout.nodes == 1) {
        // P0 traces equal the single nodal state. Use fixed-size vector views.
        for (int face = 0; face < p_.layout.cells; ++face) {
            int next = (face + 1) % p_.layout.cells;
            Eigen::Map<const FluxJacobian> left(face_left_.col(face).data());
            Eigen::Map<const FluxJacobian> right(face_right_.col(face).data());
            State flux = left * direction.segment<NVAR>(face * NVAR) +
                         right * direction.segment<NVAR>(next * NVAR);
            result.segment<NVAR>(face * NVAR) += flux;
            result.segment<NVAR>(next * NVAR) -= flux;
        }
        return;
    }
    for (int face = 0; face < p_.layout.cells; ++face) {
        int next = (face + 1) % p_.layout.cells;
        Eigen::Map<const FluxJacobian> left(face_left_.col(face).data());
        Eigen::Map<const FluxJacobian> right(face_right_.col(face).data());
        State flux = left * p_.sample(direction, face, p_.basis.right) +
                     right * p_.sample(direction, next, p_.basis.left);
        p_.layout.cell(result, face).noalias() += p_.basis.right * flux.transpose();
        p_.layout.cell(result, next).noalias() -= p_.basis.left * flux.transpose();
    }
}
void DGOperator::block(int row, int col, BlockRef result) const {
    require(ready_, "Refresh the DG Jacobian before assembling a block");
    if (row == col)
        volume_.block(row, result);
    else
        result.setZero();
    for (int side = 0; side < NFACES; ++side) {
        int face = side == 0 ? (row + p_.layout.cells - 1) % p_.layout.cells : row;
        int next = (face + 1) % p_.layout.cells;
        if (col != face && col != next)
            continue;
        const Vec &test = side == 0 ? p_.basis.left : p_.basis.right;
        double sign = side == 0 ? -1.0 : 1.0;
        // With one periodic cell both traces belong to the same block.
        for (int trace_side = 0; trace_side < NFACES; ++trace_side) {
            if (col != (trace_side == 0 ? face : next))
                continue;
            const Vec &trial = trace_side == 0 ? p_.basis.right : p_.basis.left;
            const Mat &storage = trace_side == 0 ? face_left_ : face_right_;
            Eigen::Map<const FluxJacobian> derivative(storage.col(face).data());
            for (int w = 0; w < NVAR; ++w)
                for (int j = 0; j < p_.layout.nodes; ++j)
                    for (int v = 0; v < NVAR; ++v)
                        result.col(p_.layout.local(w, j)).segment(v * p_.layout.nodes,
                                                                  p_.layout.nodes).noalias() +=
                            (sign * trial[j] * derivative(v, w)) * test;
        }
    }
}
} // namespace euler

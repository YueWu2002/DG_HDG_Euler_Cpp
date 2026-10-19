#include "euler/face_jacobian.hpp"
namespace euler {
FaceJacobian::FaceJacobian(const Problem &p)
    : p_(p), physical_left_(p.basis.nodes, NVAR * NVAR * p.face_count()),
      physical_right_(physical_left_.rows(), physical_left_.cols()),
      hll_left_(physical_left_.rows(), physical_left_.cols()),
      hll_right_(physical_left_.rows(), physical_left_.cols()),
      tau_left_(p.basis.nodes, p.face_count()), tau_right_(tau_left_.rows(), tau_left_.cols()),
      weights_(tau_left_.rows(), tau_left_.cols()), diagonal_(p.trace_size()),
      left_(p.basis.nodes, NVAR), right_(p.basis.nodes, NVAR) {
    for (int face = 0; face < p.face_count(); ++face)
        for (int q = 0; q < p.basis.nodes; ++q) weights_(q, face) = p.face_weight(face, q);
}
void FaceJacobian::refresh(const Vec &state) {
    for (int axis = 0; axis < NDIM; ++axis)
        for (int cell = 0; cell < p_.layout.cells; ++cell) {
            int face = axis * p_.layout.cells + cell, side = 2 * axis + 1;
            int next = p_.neighbor(cell, side);
            p_.face_values(p_.layout.cell(state, cell), side, left_);
            p_.face_values(p_.layout.cell(state, next), side - 1, right_);
            for (int q = 0; q < p_.basis.nodes; ++q) {
                FluxJacobian l = p_.example.gas.jacobian(left_.row(q).transpose(), axis);
                FluxJacobian r = p_.example.gas.jacobian(right_.row(q).transpose(), axis);
                auto s = p_.speeds[face * p_.basis.nodes + q];
                tau_left_(q, face) = s.a;
                tau_right_(q, face) = s.b;
                FluxJacobian gl = (weights_(q, face) * s.b / (s.a + s.b)) * (l + s.a * FluxJacobian::Identity());
                FluxJacobian gr = (weights_(q, face) * s.a / (s.a + s.b)) * (r - s.b * FluxJacobian::Identity());
                for (int w = 0; w < NVAR; ++w)
                    for (int v = 0; v < NVAR; ++v) {
                        int column = (face * NVAR + w) * NVAR + v;
                        physical_left_(q, column) = l(v, w);
                        physical_right_(q, column) = r(v, w);
                        hll_left_(q, column) = gl(v, w);
                        hll_right_(q, column) = gr(v, w);
                    }
                for (int v = 0; v < NVAR; ++v)
                    diagonal_[face * p_.face_width() + v * p_.basis.nodes + q] = -weights_(q, face) * (s.a + s.b);
            }
        }
}
void FaceJacobian::hll_action(int face, ConstBlockRef l, ConstBlockRef r, BlockRef flux) const {
    for (int v = 0; v < NVAR; ++v) {
        flux.col(v).setZero();
        for (int w = 0; w < NVAR; ++w) {
            int column = (face * NVAR + w) * NVAR + v;
            flux.col(v).array() += hll_left_.col(column).array() * l.col(w).array() +
                                   hll_right_.col(column).array() * r.col(w).array();
        }
    }
}
FluxJacobian FaceJacobian::hll(int face, int q, bool left) const {
    const Mat &storage = left ? hll_left_ : hll_right_;
    FluxJacobian result;
    for (int w = 0; w < NVAR; ++w)
        for (int v = 0; v < NVAR; ++v) result(v, w) = storage(q, (face * NVAR + w) * NVAR + v);
    return result;
}
double FaceJacobian::tau(int face, int q, bool left) const {
    return left ? tau_left_(q, face) : tau_right_(q, face);
}
FluxJacobian FaceJacobian::hdg(int face, int q, bool left) const {
    const Mat &storage = left ? physical_left_ : physical_right_;
    FluxJacobian result;
    for (int w = 0; w < NVAR; ++w)
        for (int v = 0; v < NVAR; ++v)
            result(v, w) = (left ? 1 : -1) * storage(q, (face * NVAR + w) * NVAR + v);
    result.diagonal().array() += tau(face, q, left);
    result *= weights_(q, face);
    return result;
}
void FaceJacobian::cell_c(int cell, Eigen::Ref<const Vec> direction, Eigen::Ref<Vec> result) const {
    Eigen::Map<const Mat> nodal(direction.data(), p_.layout.plane(), NVAR);
    auto fs = p_.faces(cell);
    for (int side = 0; side < NFACES; ++side) {
        int face = fs[side];
        bool left = side % 2 != 0;
        const Mat &physical = left ? physical_left_ : physical_right_;
        const Mat &stabilization = left ? tau_left_ : tau_right_;
        p_.face_values(nodal, side, left_);
        Eigen::Map<Mat> output(result.data() + side * p_.face_width(), p_.basis.nodes, NVAR);
        for (int v = 0; v < NVAR; ++v) {
            output.col(v).array() = stabilization.col(face).array() * left_.col(v).array();
            for (int w = 0; w < NVAR; ++w)
                output.col(v).array() += (left ? 1 : -1) * physical.col((face * NVAR + w) * NVAR + v).array() * left_.col(w).array();
            output.col(v).array() *= weights_.col(face).array();
        }
    }
}
} // namespace euler

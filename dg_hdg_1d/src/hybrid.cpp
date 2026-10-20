#include "euler/hybrid.hpp"
namespace euler {
HybridOperator::HybridOperator(const Problem &p)
    : p_(p), volume_(p), volume_residual_(p.layout.size()), skeleton_residual_(p.layout.cells * NVAR),
      rhs_(p.layout.cells * NVAR), local_(p.layout.block_size()), faces_(NFACES * NVAR),
      face_result_(NFACES * NVAR), schur_error_(rhs_.size()) {
    cells_.reserve(p.layout.cells);
    for (int c = 0; c < p.layout.cells; ++c)
        cells_.emplace_back(p.layout.block_size());
}
int HybridOperator::face_index(int cell, int side) const {
    return side == 0 ? (cell + p_.layout.cells - 1) % p_.layout.cells : cell;
}
void HybridOperator::gather(const Vec &x, int cell) const {
    for (int side = 0; side < NFACES; ++side)
        faces_.segment<NVAR>(side * NVAR) = x.segment<NVAR>(face_index(cell, side) * NVAR);
}
void HybridOperator::residual(const Vec &u, const Vec &lambda, Vec &r, Vec &g) const {
    p_.volume_residual(u, r);
    for (int face = 0; face < p_.layout.cells; ++face) {
        int next = (face + 1) % p_.layout.cells;
        auto s = p_.speeds[face];
        State left = p_.sample(u, face, p_.basis.right), right = p_.sample(u, next, p_.basis.left);
        State z = lambda.segment<NVAR>(face * NVAR);
        State fl = p_.example.gas.flux(left) + s.a * (left - z);
        State fr = p_.example.gas.flux(right) + s.b * (z - right);
        g.segment<NVAR>(face * NVAR) = fl - fr;
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < p_.layout.nodes; ++i) {
                r[p_.layout.index(face, v, i)] += p_.basis.right[i] * fl[v];
                r[p_.layout.index(next, v, i)] -= p_.basis.left[i] * fr[v];
            }
    }
}
void HybridOperator::refresh(const Vec &u, const Vec &lambda) {
    residual(u, lambda, volume_residual_, skeleton_residual_);
    condense(u);
}
void HybridOperator::refresh_projected(const Vec &u, const Vec &r) {
    volume_residual_ = r;
    skeleton_residual_.setZero();
    condense(u);
}
void HybridOperator::condense(const Vec &u) {
    volume_.refresh(u);
    rhs_ = -skeleton_residual_;
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        auto &data = cells_[cell];
        volume_.block(cell, data.A);
        data.B.setZero();
        data.C.setZero();
        for (int side = 0; side < NFACES; ++side) {
            auto speed = p_.speeds[face_index(cell, side)];
            const Vec &shape = side == 0 ? p_.basis.left : p_.basis.right;
            double sign = side == 0 ? -1.0 : 1.0;
            double stabilization = side == 0 ? speed.b : speed.a;
            FluxJacobian derivative = sign * p_.example.gas.jacobian(p_.sample(u, cell, shape)) +
                                      stabilization * FluxJacobian::Identity();
            for (int v = 0; v < NVAR; ++v)
                for (int i = 0; i < p_.layout.nodes; ++i) {
                    int row = p_.layout.local(v, i);
                    data.B(row, side * NVAR + v) = -shape[i] * stabilization;
                    for (int w = 0; w < NVAR; ++w) {
                        data.C(side * NVAR + w, row) = shape[i] * derivative(w, v);
                        for (int j = 0; j < p_.layout.nodes; ++j)
                            data.A(row, p_.layout.local(w, j)) +=
                                shape[i] * shape[j] * derivative(v, w);
                    }
                }
        }
        data.lu.compute(data.A);
        for (int j = 0; j < NFACES * NVAR; ++j)
            data.lu.solve(data.B.col(j), data.X.col(j));
        data.lu.solve(
            volume_residual_.segment(cell * p_.layout.block_size(), p_.layout.block_size()),
            data.y);
        for (int j = 0; j < NFACES * NVAR; ++j)
            data.S.col(j).noalias() = -data.C * data.X.col(j);
        face_result_.noalias() = data.C * data.y;
        for (int side = 0; side < NFACES; ++side)
            rhs_.segment<NVAR>(face_index(cell, side) * NVAR) +=
                face_result_.segment<NVAR>(side * NVAR);
    }
}
void HybridOperator::apply(const Vec &x, Vec &y) const {
    for (int face = 0; face < p_.layout.cells; ++face)
        y.segment<NVAR>(face * NVAR) =
            -(p_.speeds[face].a + p_.speeds[face].b) * x.segment<NVAR>(face * NVAR);
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        gather(x, cell);
        // In 1D the condensed cell block is always 6x6, independent of k.
        // Fixed Eigen maps keep this tiny multiply inlined instead of making
        // one dynamic BLAS GEMV call per cell and Krylov iteration.
        using FaceVector = Eigen::Matrix<double, NFACES * NVAR, 1>;
        using FaceMatrix = Eigen::Matrix<double, NFACES * NVAR, NFACES * NVAR>;
        Eigen::Map<FaceVector>(face_result_.data()).noalias() =
            Eigen::Map<const FaceMatrix>(cells_[cell].S.data()) *
            Eigen::Map<const FaceVector>(faces_.data());
        for (int side = 0; side < NFACES; ++side)
            y.segment<NVAR>(face_index(cell, side) * NVAR) +=
                face_result_.segment<NVAR>(side * NVAR);
    }
}
void HybridOperator::block(int row, int col, BlockRef result) const {
    result.setZero();
    if (row == col)
        result.diagonal().setConstant(-(p_.speeds[row].a + p_.speeds[row].b));
    // A face belongs to the cell on its left and the cell on its right.
    for (int adjacent = 0; adjacent < NFACES; ++adjacent) {
        int cell = (row + adjacent) % p_.layout.cells;
        const auto &data = cells_[cell];
        int row_side = adjacent == 0 ? 1 : 0;
        for (int col_side = 0; col_side < NFACES; ++col_side)
            if (face_index(cell, col_side) == col)
                result += data.S.block<NVAR, NVAR>(row_side * NVAR, col_side * NVAR);
    }
}
void HybridOperator::recovered_residual(const Vec &dl, Vec &result) const {
    apply(dl, schur_error_);
    schur_error_ -= rhs_;
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        gather(schur_error_, cell);
        for (int side = 0; side < NFACES; ++side) {
            const auto s = p_.speeds[face_index(cell, side)];
            faces_.segment<NVAR>(side * NVAR) /= -(s.a + s.b);
        }
        result.segment(cell * p_.layout.block_size(), p_.layout.block_size()).noalias() =
            -cells_[cell].B * faces_;
    }
}
void HybridOperator::recover(const Vec &dl, Vec &du) const {
    for (int cell = 0; cell < p_.layout.cells; ++cell) {
        gather(dl, cell);
        local_.noalias() = cells_[cell].X * faces_;
        du.segment(cell * p_.layout.block_size(), p_.layout.block_size()) =
            -cells_[cell].y - local_;
    }
}
} // namespace euler

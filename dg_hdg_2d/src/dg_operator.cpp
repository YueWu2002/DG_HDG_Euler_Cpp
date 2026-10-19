#include "euler/dg_operator.hpp"
namespace euler {
DGOperator::DGOperator(const Problem &p)
    : p_(p), volume_(p), faces_(p), left_(p.basis.nodes, NVAR),
      right_(p.basis.nodes, NVAR), flux_(p.basis.nodes, NVAR) {}
void DGOperator::refresh(const Vec &state) {
    volume_.refresh(state);
    faces_.refresh(state);
}
void DGOperator::apply(const Vec &direction, Vec &result) const {
    volume_.apply(direction, result);
    for (int axis = 0; axis < NDIM; ++axis)
        for (int cell = 0; cell < p_.layout.cells; ++cell) {
            int side = 2 * axis + 1, next = p_.neighbor(cell, side), face = axis * p_.layout.cells + cell;
            p_.face_values(p_.layout.cell(direction, cell), side, left_);
            p_.face_values(p_.layout.cell(direction, next), side - 1, right_);
            faces_.hll_action(face, left_, right_, flux_);
            p_.lift(side, flux_, p_.layout.cell(result, cell));
            p_.lift(side - 1, flux_, p_.layout.cell(result, next), -1);
        }
}
void DGOperator::block(int row, int col, BlockRef result) const {
    if (row == col) volume_.block(row, result); else result.setZero();
    auto fs = p_.faces(row);
    for (int side = 0; side < NFACES; ++side) {
        if (col != row && col != p_.neighbor(row, side)) continue;
        bool left = col == (fs[side] % p_.layout.cells);
        int trial_side = 2 * (side / 2) + (left ? 1 : 0);
        const Vec &test = side % 2 ? p_.basis.right : p_.basis.left;
        const Vec &trial = left ? p_.basis.right : p_.basis.left;
        double sign = side % 2 ? 1 : -1;
        for (int q = 0; q < p_.basis.nodes; ++q) {
            FluxJacobian derivative = faces_.hll(fs[side], q, left);
            for (int w = 0; w < NVAR; ++w)
                for (int j = 0; j < p_.basis.nodes; ++j) {
                    auto column = result.col(p_.layout.local(w, p_.face_node(trial_side, q, j)));
                    for (int v = 0; v < NVAR; ++v) {
                        double scale = sign * trial[j] * derivative(v, w);
                        if (side < 2)
                            column.segment(p_.layout.local(v, q * p_.basis.nodes), p_.basis.nodes).noalias() += scale * test;
                        else
                            for (int i = 0; i < p_.basis.nodes; ++i)
                                column[p_.layout.local(v, p_.face_node(side, q, i))] += scale * test[i];
                    }
                }
        }
    }
}
} // namespace euler

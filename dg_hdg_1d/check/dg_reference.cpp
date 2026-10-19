#include "dg_reference.hpp"
namespace euler::check::reference {
namespace {
State quadrature_state(const Problem &p, const Vec &u, int cell, int q) {
    State result = State::Zero();
    for (int v = 0; v < NVAR; ++v)
        for (int i = 0; i < p.layout.nodes; ++i)
            result[v] += p.basis.values(q, i) * u[p.layout.index(cell, v, i)];
    return result;
}
} // namespace
void volume_block(const Problem &p, const Vec &u, int e, BlockRef a) {
    a.setZero();
    a.diagonal() = p.mass.segment(e * p.layout.block_size(), p.layout.block_size()) / p.step_size;
    for (int q = 0; q < p.basis.points; ++q) {
        FluxJacobian f = p.example.gas.jacobian(quadrature_state(p, u, e, q));
        for (int v = 0; v < NVAR; ++v)
            for (int w = 0; w < NVAR; ++w)
                for (int i = 0; i < p.layout.nodes; ++i)
                    for (int j = 0; j < p.layout.nodes; ++j)
                        a(p.layout.local(v, i), p.layout.local(w, j)) -=
                            p.basis.weighted_derivatives(q, i) * p.basis.values(q, j) * f(v, w);
    }
}
void dg_action(const Problem &p, const Vec &u, const Vec &x, Vec &y) {
    y.array() = p.mass.array() * x.array() / p.step_size;
    for (int e = 0; e < p.layout.cells; ++e)
        for (int q = 0; q < p.basis.points; ++q) {
            State f = p.example.gas.jacobian(quadrature_state(p, u, e, q)) *
                      quadrature_state(p, x, e, q);
            for (int v = 0; v < NVAR; ++v)
                for (int i = 0; i < p.layout.nodes; ++i)
                    y[p.layout.index(e, v, i)] -= p.basis.weighted_derivatives(q, i) * f[v];
        }
    for (int e = 0; e < p.layout.cells; ++e) {
        int next = (e + 1) % p.layout.cells;
        auto s = p.speeds[e];
        State xl = p.sample(x, e, p.basis.right), xr = p.sample(x, next, p.basis.left);
        State f =
            (s.b * p.example.gas.jacobian(p.sample(u, e, p.basis.right)) * xl +
             s.a * p.example.gas.jacobian(p.sample(u, next, p.basis.left)) * xr + s.a * s.b * (xl - xr)) /
            (s.a + s.b);
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < p.layout.nodes; ++i) {
                y[p.layout.index(e, v, i)] += p.basis.right[i] * f[v];
                y[p.layout.index(next, v, i)] -= p.basis.left[i] * f[v];
            }
    }
}
} // namespace euler::check::reference

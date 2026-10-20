#include "euler/discretization.hpp"
#include <algorithm>
namespace euler {
Problem::Problem(const Config &c)
    : config(c), basis(c.degree, c.quadrature), layout{c.cells, basis.nodes}, h(2.0 / c.cells),
      old(layout.size()), mass(layout.size()), inverse_mass(layout.size()), speeds(c.cells),
      derivative_transpose_(basis.weighted_derivatives.transpose()), flux_workspace_(basis.points, NVAR) {
    validate(c);
    for (int cell = 0; cell < layout.cells; ++cell)
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i)
                mass[layout.index(cell, v, i)] = .5 * h * basis.w[i];
    inverse_mass = mass.cwiseInverse();
}
State Problem::sample(const Vec &u, int cell, const Vec &shape) const {
    State result = State::Zero();
    for (int v = 0; v < NVAR; ++v)
        for (int i = 0; i < layout.nodes; ++i)
            result[v] += shape[i] * u[layout.index(cell, v, i)];
    return result;
}
State Problem::quadrature_state(const Vec &u, int cell, int q) const {
    if (basis.points == layout.nodes)
        return layout.cell(u, cell).row(q).transpose();
    State result = State::Zero();
    for (int v = 0; v < NVAR; ++v)
        for (int i = 0; i < layout.nodes; ++i)
            result[v] += basis.values(q, i) * u[layout.index(cell, v, i)];
    return result;
}
void Problem::initial(Vec &u) const {
    Mat samples(basis.error_x.size(), NVAR);
    Mat nodal(layout.nodes, NVAR);
    for (int cell = 0; cell < layout.cells; ++cell) {
        for (int q = 0; q < basis.error_x.size(); ++q)
            samples.row(q) =
                example.exact(-1 + h * (cell + .5) + .5 * h * basis.error_x[q], 0).transpose();
        nodal.noalias() = basis.projection * samples;
        State mean = .5 * nodal.transpose() * basis.w;
        require(mean[Density] > 0 && example.gas.pressure(mean) > 0, "Nonphysical cell average");
        // Convex scaling preserves the cell integral. Include all points used by
        // the residual, boundary flux, and error integration.
        Mat checks(layout.nodes + basis.points + basis.error_x.size() + NFACES, NVAR);
        checks.topRows(layout.nodes) = nodal;
        checks.middleRows(layout.nodes, basis.points).noalias() = basis.values * nodal;
        checks.middleRows(layout.nodes + basis.points, basis.error_x.size()).noalias() =
            basis.error_values * nodal;
        checks.row(checks.rows() - NFACES) = basis.left.transpose() * nodal;
        checks.row(checks.rows() - 1) = basis.right.transpose() * nodal;
        double theta = 1;
        for (int q = 0; q < checks.rows(); ++q) {
            State difference = checks.row(q).transpose() - mean;
            auto physical = [&](double factor) {
                State z = mean + factor * difference;
                return z[Density] > 1e-12 && example.gas.pressure(z) > 1e-12;
            };
            if (!physical(theta)) {
                double low = 0, high = theta;
                for (int k = 0; k < 60; ++k) {
                    double mid = .5 * (low + high);
                    if (physical(mid))
                        low = mid;
                    else
                        high = mid;
                }
                theta = .99 * low;
            }
        }
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i)
                u[layout.index(cell, v, i)] = mean[v] + theta * (nodal(i, v) - mean[v]);
    }
    require(admissible(u), "Invalid projected initial state");
}
double Problem::freeze(const Vec &u) {
    old = u;
    double rate = 0;
    for (int e = 0; e < layout.cells; ++e)
        speeds[e] =
            einfeldt(sample(u, e, basis.right), sample(u, (e + 1) % layout.cells, basis.left),
                     example.gas, config.safety);
    for (int e = 0; e < layout.cells; ++e) {
        const auto &l = speeds[(e + layout.cells - 1) % layout.cells];
        const auto &r = speeds[e];
        double alpha = std::max({l.a, l.b, r.a, r.b});
        for (int q = 0; q < basis.points; ++q)
            alpha = std::max(alpha, example.gas.spectral_radius(quadrature_state(u, e, q)));
        rate = std::max(rate, alpha / h);
    }
    require(std::isfinite(rate) && rate > 0, "Invalid CFL rate");
    return rate;
}
void Problem::trace(const Vec &u, Vec &lambda) const {
    for (int face = 0; face < layout.cells; ++face) {
        State l = sample(u, face, basis.right),
              r = sample(u, (face + 1) % layout.cells, basis.left);
        auto s = speeds[face];
        lambda.segment<NVAR>(face * NVAR) =
            (s.a * l + s.b * r + example.gas.flux(l) - example.gas.flux(r)) / (s.a + s.b);
    }
}
void Problem::volume_residual(const Vec &u, Vec &r) const {
    r.array() = mass.array() * (u - old).array() / step_size;
    if (basis.nodes == 1) return;
    for (int e = 0; e < layout.cells; ++e) {
        for (int q = 0; q < basis.points; ++q)
            flux_workspace_.row(q) = example.gas.flux(quadrature_state(u, e, q)).transpose();
        auto cell = layout.cell(r, e);
        for (int v = 0; v < NVAR; ++v)
            cell.col(v).noalias() -= derivative_transpose_ * flux_workspace_.col(v);
    }
}
void Problem::dg_residual(const Vec &u, Vec &r) const {
    volume_residual(u, r);
    for (int e = 0; e < layout.cells; ++e) {
        int next = (e + 1) % layout.cells;
        auto s = speeds[e];
        State l = sample(u, e, basis.right), rr = sample(u, next, basis.left);
        State f = (s.b * example.gas.flux(l) + s.a * example.gas.flux(rr) + s.a * s.b * (l - rr)) /
                  (s.a + s.b);
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i) {
                r[layout.index(e, v, i)] += basis.right[i] * f[v];
                r[layout.index(next, v, i)] -= basis.left[i] * f[v];
            }
    }
}
bool Problem::admissible(const Vec &u) const {
    if (!u.allFinite())
        return false;
    auto ok = [&](const State &z) { return z[Density] > 1e-13 && example.gas.pressure(z) > 1e-13; };
    for (int e = 0; e < layout.cells; ++e) {
        auto cell = layout.cell(u, e);
        for (int i = 0; i < layout.nodes; ++i)
            if (!ok(cell.row(i).transpose()))
                return false;
        for (int q = 0; q < basis.points; ++q)
            if (!ok(quadrature_state(u, e, q)))
                return false;
        if (!ok(sample(u, e, basis.left)) || !ok(sample(u, e, basis.right)))
            return false;
    }
    return true;
}
double Problem::residual_norm(const Vec &r) const {
    return std::sqrt(.5 * (r.array().square() * inverse_mass.array()).sum());
}
State Problem::integral(const Vec &u) const {
    State total = State::Zero();
    for (int e = 0; e < layout.cells; ++e)
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i)
                total[v] += .5 * h * basis.w[i] * u[layout.index(e, v, i)];
    return total;
}
State Problem::error(const Vec &u, double t) const {
    State squared = State::Zero();
    for (int e = 0; e < layout.cells; ++e)
        for (int q = 0; q < basis.error_x.size(); ++q) {
            State z = State::Zero();
            for (int v = 0; v < NVAR; ++v)
                for (int i = 0; i < layout.nodes; ++i)
                    z[v] += basis.error_values(q, i) * u[layout.index(e, v, i)];
            State difference = z - example.exact(-1 + h * (e + .5) + .5 * h * basis.error_x[q], t);
            squared.array() += .5 * h * basis.error_w[q] * difference.array().square();
        }
    return squared.cwiseSqrt();
}
} // namespace euler

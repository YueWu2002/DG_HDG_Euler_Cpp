#include "euler/discretization.hpp"
#include "euler/tensor_kernels.hpp"
namespace euler {
void Problem::initial(Vec &u) const {
    int n = basis.nodes, Q = static_cast<int>(basis.error_x.size());
    NodalBasis error_basis(config.degree, Q);
    TensorProduct error_tensor(error_basis, hx, hy);
    Mat samples(Q * Q, NVAR), scratch(n, Q), points(basis.points * basis.points, NVAR);
    Mat checks(n * n + points.rows() + Q * Q + NFACES * n, NVAR), values(Q * Q, NVAR);
    for (int cell = 0; cell < layout.cells; ++cell) {
        for (int j = 0; j < Q; ++j)
            for (int i = 0; i < Q; ++i)
                samples.row(i + Q * j) = example.exact(
                    -.5 * example.lx + hx * (cell % layout.nx + .5) + .5 * hx * basis.error_x[i],
                    -.5 * example.ly + hy * (cell / layout.nx + .5) + .5 * hy * basis.error_x[j], 0).transpose();
        for (int v = 0; v < NVAR; ++v) {
            Eigen::Map<const Mat> sample(samples.col(v).data(), Q, Q);
            multiply(basis.projection, sample, scratch);
            multiply_transpose_right(scratch, basis.projection, layout.sheet(u, cell, v));
        }
        auto nodal = layout.cell(u, cell);
        State mean = State::Zero();
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                mean += .25 * basis.w[i] * basis.w[j] * nodal.row(i + n * j).transpose();
        require(mean[Density] > 0 && example.gas.pressure(mean) > 0, "Nonphysical cell average");
        checks.topRows(n * n) = nodal;
        tensor_.interpolate(nodal, points);
        checks.middleRows(n * n, points.rows()) = points;
        error_tensor.interpolate(nodal, values);
        checks.middleRows(n * n + points.rows(), Q * Q) = values;
        for (int side = 0; side < NFACES; ++side)
            face_values(nodal, side, checks.middleRows(n * n + points.rows() + Q * Q + side * n, n));
        double theta = 1;
        for (int q = 0; q < checks.rows(); ++q) {
            State difference = checks.row(q).transpose() - mean;
            auto physical = [&](double factor) {
                State z = mean + factor * difference;
                return z[Density] > 1e-12 && example.gas.pressure(z) > 1e-12;
            };
            if (!physical(theta)) {
                double low = 0, high = theta;
                for (int it = 0; it < 60; ++it) {
                    double middle = .5 * (low + high);
                    if (physical(middle)) low = middle; else high = middle;
                }
                theta = .99 * low;
            }
        }
        for (int v = 0; v < NVAR; ++v)
            nodal.col(v).array() = mean[v] + theta * (nodal.col(v).array() - mean[v]);
    }
    require(admissible(u), "Invalid projected initial state");
}
State Problem::error(const Vec &u, double time) const {
    int Q = static_cast<int>(basis.error_x.size());
    NodalBasis error_basis(config.degree, Q);
    TensorProduct error_tensor(error_basis, hx, hy);
    Mat values(Q * Q, NVAR);
    State result = State::Zero();
    for (int cell = 0; cell < layout.cells; ++cell) {
        error_tensor.interpolate(layout.cell(u, cell), values);
        for (int j = 0; j < Q; ++j)
            for (int i = 0; i < Q; ++i) {
                State difference = values.row(i + Q * j).transpose() - example.exact(
                    -.5 * example.lx + hx * (cell % layout.nx + .5) + .5 * hx * basis.error_x[i],
                    -.5 * example.ly + hy * (cell / layout.nx + .5) + .5 * hy * basis.error_x[j], time);
                result.array() += .25 * hx * hy * basis.error_w[i] * basis.error_w[j] *
                                  difference.array().square();
            }
    }
    return result.array().sqrt();
}
} // namespace euler

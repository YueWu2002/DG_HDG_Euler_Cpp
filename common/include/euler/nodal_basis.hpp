#pragma once
#include "types.hpp"
#include <fastgl.h>
#include <algorithm>
#include <utility>
#include <tuple>
namespace euler {
// Setup-only, header-only reference element. No wrapper basis.cpp.
struct NodalBasis {
    int nodes, points;
    Vec x, w, qx, qw, left, right, error_x, error_w;
    Mat values, weighted_derivatives, error_values, projection;
    static std::pair<Vec, Vec> gauss(int n) {
        require(n > 0, "Gauss rule requires a positive point count");
        Vec nodes(n), weights(n);
        for (int k = 1; k <= (n + 1) / 2; ++k) {
            auto pair = fastgl::GLPair(n, k);
            nodes[k - 1] = -pair.x();
            nodes[n - k] = pair.x();
            weights[k - 1] = weights[n - k] = pair.weight;
        }
        if (n % 2)
            nodes[n / 2] = 0;
        return {std::move(nodes), std::move(weights)};
    }
    Vec evaluate(double z) const {
        Vec v = Vec::Ones(nodes);
        for (int i = 0; i < nodes; ++i)
            for (int j = 0; j < nodes; ++j)
                if (i != j)
                    v[i] *= (z - x[j]) / (x[i] - x[j]);
        return v;
    }
    explicit NodalBasis(int degree, int quadrature = 0)
        : nodes(degree + 1), points(quadrature ? quadrature : nodes) {
        require(nodes > 0 && points >= nodes, "Require quadrature >= degree+1");
        std::tie(x, w) = gauss(nodes);
        std::tie(qx, qw) = gauss(points);
        Vec barycentric = Vec::Ones(nodes);
        for (int i = 0; i < nodes; ++i)
            for (int j = 0; j < nodes; ++j)
                if (i != j)
                    barycentric[i] /= x[i] - x[j];
        Mat derivative = Mat::Zero(nodes, nodes);
        for (int i = 0; i < nodes; ++i) {
            for (int j = 0; j < nodes; ++j)
                if (i != j)
                    derivative(i, j) = barycentric[j] / (barycentric[i] * (x[i] - x[j]));
            derivative(i, i) = -derivative.row(i).sum();
        }
        values.resize(points, nodes);
        for (int q = 0; q < points; ++q)
            values.row(q) = evaluate(qx[q]).transpose();
        if (points == nodes)
            values.setIdentity();
        weighted_derivatives = qw.asDiagonal() * values * derivative;
        left = evaluate(-1);
        right = evaluate(1);
        std::tie(error_x, error_w) = gauss(std::max(20, 2 * nodes + 3));
        error_values.resize(error_x.size(), nodes);
        for (int q = 0; q < error_x.size(); ++q)
            error_values.row(q) = evaluate(error_x[q]).transpose();
        projection =
            w.cwiseInverse().asDiagonal() * error_values.transpose() * error_w.asDiagonal();
    }
};
} // namespace euler

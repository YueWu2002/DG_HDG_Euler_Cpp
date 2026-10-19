#pragma once
#include "physics.hpp"
namespace euler {
// This example owns its gamma; it is not a CLI or per-step output field.
struct SimpleWave {
    Gas gas{3.0};
    State exact(double x, double time) const;
    double shock_time() const { return 4 / ((gas.gamma() + 1) * std::acos(-1.0)); }
};
} // namespace euler

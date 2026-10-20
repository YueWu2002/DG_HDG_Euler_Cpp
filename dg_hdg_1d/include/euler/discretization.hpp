#pragma once
#include "config.hpp"
#include "euler/nodal_basis.hpp"
#include "simple_wave.hpp"
#include <vector>
namespace euler {
class Problem {
  public:
    Config config;
    SimpleWave example;
    NodalBasis basis;
    Layout layout;
    double h, step_size = 0;
    Vec old, mass, inverse_mass;
    std::vector<Speed> speeds;
    explicit Problem(const Config &config);
    void initial(Vec &u) const;
    State sample(const Vec &u, int cell, const Vec &shape) const;
    State quadrature_state(const Vec &u, int cell, int q) const;
    double freeze(const Vec &u); // freeze U^n and speeds; return CFL spectral rate
    void trace(const Vec &u, Vec &lambda) const;
    void volume_residual(const Vec &u, Vec &residual) const;
    void dg_residual(const Vec &u, Vec &residual) const;
    bool admissible(const Vec &u) const;
    double residual_norm(const Vec &r) const;
    State integral(const Vec &u) const;
    State error(const Vec &u, double time) const;

  private:
    Mat derivative_transpose_;
    mutable Mat flux_workspace_;
};
} // namespace euler

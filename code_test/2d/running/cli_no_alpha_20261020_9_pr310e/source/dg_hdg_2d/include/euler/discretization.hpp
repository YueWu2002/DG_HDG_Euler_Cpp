#pragma once
#include "config.hpp"
#include "example.hpp"
#include "tensor_product.hpp"
#include "euler/linalg/block_graph.hpp"
namespace euler {
class Problem {
  public:
    Config config;
    NodalBasis basis;
    Layout layout;
    Example example;
    double hx, hy, step_size = 0;
    Vec old, mass, inverse_mass, trace_mass, inverse_trace_mass;
    std::vector<Speed> speeds;
    explicit Problem(const Config &config);
    int neighbor(int cell, int side) const;
    std::array<int, NFACES> faces(int cell) const;
    int face_width() const { return basis.nodes * NVAR; }
    int face_count() const { return NDIM * layout.cells; }
    int trace_size() const { return face_count() * face_width(); }
    int local_trace_size() const { return NFACES * face_width(); }
    double face_weight(int face, int q) const;
    int face_node(int side, int tangent, int normal) const;
    State sample(const Vec &u, int cell, int side, int q) const;
    void face_values(ConstBlockRef nodal, int side, BlockRef values) const;
    void lift(int side, ConstBlockRef values, BlockRef nodal, double scale = 1) const;
    void initial(Vec &state) const;
    double freeze(const Vec &state);
    void trace(const Vec &state, Vec &lambda) const;
    void volume_residual(const Vec &state, Vec &result) const;
    void dg_residual(const Vec &state, Vec &result) const;
    void hybrid_residual(const Vec &state, const Vec &lambda, Vec &volume, Vec &skeleton) const;
    bool admissible(const Vec &state) const;
    double residual_norm(const Vec &residual) const;
    State integral(const Vec &state) const;
    State error(const Vec &state, double time) const;
    BlockGraph trace_graph() const;

  private:
    TensorProduct tensor_;
    mutable Mat points_, flux_x_, flux_y_, left_, right_, face_flux_, other_flux_;
    void evaluate_volume_flux(ConstBlockRef states) const;
    bool physical(ConstBlockRef states) const;
};
} // namespace euler

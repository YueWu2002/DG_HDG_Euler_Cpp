#pragma once
#include "euler/types.hpp"
namespace euler {
// Non-owning, allocation-free callback; no std::function capture allocations.
struct Action {
    const void *context;
    void (*call)(const void *, const Vec &, Vec &);
    void operator()(const Vec &x, Vec &y) const { call(context, x, y); }
    template <class Operator> static Action bind(const Operator &op) {
        return {&op, [](const void *p, const Vec &x, Vec &y) {
                    static_cast<const Operator *>(p)->apply(x, y);
                }};
    }
};
struct LinearResult {
    int iterations;
    double relative_residual;
};
// Restarted right-preconditioned GMRES with twice classical Gram-Schmidt.
// A and M^-1 are operator actions. Eigen owns all numerical workspaces.
class GMRES {
  public:
    GMRES(int size, int restart, const Vec &inverse_mass);
    LinearResult solve(Action matrix, Action preconditioner, const Vec &rhs, Vec &x,
                       int max_iterations, double tolerance);

    const Vec &residual() const { return residual_; }

  private:
    int restart_;
    Mat basis_, hessenberg_;
    Vec cosine_, sine_, rotated_rhs_, coefficients_, work_, residual_, column_, weights_, weighted_, projection_;
    double dot(const Vec &a, const Vec &b) const;
    double norm(const Vec &a) const;
};
} // namespace euler

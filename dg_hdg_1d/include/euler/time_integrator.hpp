#pragma once
#include "hybrid.hpp"
#include "dg_operator.hpp"
#include "linalg/block_ilu.hpp"
#include "linalg/gmres.hpp"
#include "linalg/direct.hpp"
#include <optional>
namespace euler {
struct RunStatistics {
    double final_time = 0, conservation_drift = 0, seconds = 0, max_linear_residual = 0;
    double minimum_step = 0, maximum_step = 0;
    int steps = 0, newton_iterations = 0, krylov_iterations = 0;
};
// All capacity is acquired in the constructor, before advance(). The object is
// reusable for another run on the same mesh; no resize occurs during stepping.
class TimeIntegrator {
  public:
    TimeIntegrator(Problem &problem, Method method);
    RunStatistics advance(Vec &solution, bool forbid_eigen_allocations = true);

  private:
    Problem &p_;
    Method method_;
    Vec state_, trial_, residual_, trial_residual_, direction_, rhs_;
    Vec lambda_, trial_lambda_, skeleton_, trial_skeleton_, trace_direction_;
    std::optional<HybridOperator> hybrid_;
    std::optional<DGOperator> dg_;
    std::optional<BlockPreconditioner> preconditioner_;
    std::optional<GMRES> gmres_;
    std::optional<DirectSolver> direct_;
    Vec linear_weights_;
    double evaluate(const Vec &u, const Vec &lambda, Vec &r, Vec &g);
};
} // namespace euler

#pragma once
#include "hybrid.hpp"
#include "dg_operator.hpp"
#include "euler/linalg/block_jacobi.hpp"
#include "euler/linalg/gmres.hpp"
#include "euler/linalg/direct.hpp"
#include <optional>
namespace euler {
struct UpdateStatistics {
    int step, stage, iteration, krylov_iterations;
    double time, stage_time, dt, cfl, residual_physical_rms, damping;
    double linear_relative_residual, dg_linear_relative_residual, iteration_s;
};
struct RunStatistics {
    double final_time = 0, conservation_drift = 0, seconds = 0, max_linear_residual = 0;
    double max_dg_linear_residual = 0, minimum_step = 0, maximum_step = 0;
    double measured_newton_seconds = 0;
    int steps = 0, newton_iterations = 0, krylov_iterations = 0, measured_updates = 0;
    std::vector<UpdateStatistics> updates;
};
class TimeIntegrator {
  public:
    TimeIntegrator(Problem &problem, Method method);
    RunStatistics advance(Vec &solution, bool forbid_eigen_allocations = true);

  private:
    Problem &p_;
    Method method_;
    Vec state_, trial_, residual_, trial_residual_, direction_, rhs_;
    Vec lambda_, trace_direction_, equation_error_, independent_error_, step_state_, stage_base_;
    Mat slopes_;
    std::optional<HybridOperator> hybrid_;
    std::optional<DGOperator> dg_, diagnostic_;
    std::optional<BlockJacobi> preconditioner_;
    std::optional<GMRES> gmres_;
    std::optional<DirectSolver> direct_;
    Vec linear_weights_;
    double dual_norm(const Vec &x) const;
    double state_norm(const Vec &x) const;
};
} // namespace euler

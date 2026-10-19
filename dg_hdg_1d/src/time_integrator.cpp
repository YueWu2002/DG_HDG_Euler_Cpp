#include "euler/time_integrator.hpp"
#include <algorithm>
#include <chrono>
#include <limits>
namespace euler {
namespace {
// Restore the caller's allocation state even if an iteration throws.
class AllocationGuard {
    bool previous_;

  public:
    explicit AllocationGuard(bool forbid) : previous_(Eigen::internal::is_malloc_allowed()) {
        if (forbid)
            Eigen::internal::set_is_malloc_allowed(false);
    }
    ~AllocationGuard() { Eigen::internal::set_is_malloc_allowed(previous_); }
};
} // namespace
TimeIntegrator::TimeIntegrator(Problem &p, Method method)
    : p_(p), method_(method), state_(p.layout.size()), trial_(p.layout.size()),
      residual_(p.layout.size()), trial_residual_(p.layout.size()), direction_(p.layout.size()),
      rhs_(method == Method::DG ? p.layout.size() : p.layout.cells * NVAR),
      lambda_(p.layout.cells * NVAR), trial_lambda_(p.layout.cells * NVAR),
      skeleton_(p.layout.cells * NVAR), trial_skeleton_(p.layout.cells * NVAR),
      trace_direction_(p.layout.cells * NVAR),
      linear_weights_(method == Method::DG ? p.inverse_mass : Vec::Ones(rhs_.size())) {
    if (method == Method::HDG)
        hybrid_.emplace(p);
    else
        dg_.emplace(p);
    auto graph = BlockGraph::periodic(p.layout.cells);
    int width = method == Method::DG ? p.layout.block_size() : NVAR;
    if (p.config.linear == LinearSolver::Direct)
        direct_.emplace(graph, width);
    else {
        preconditioner_.emplace(graph, width, p.config.preconditioner, p.config.colamd);
        gmres_.emplace(rhs_.size(), p.config.restart, linear_weights_);
    }
}
double TimeIntegrator::evaluate(const Vec &u, const Vec &lambda, Vec &r, Vec &g) {
    if (method_ == Method::DG) {
        p_.dg_residual(u, r);
        return p_.residual_norm(r);
    }
    hybrid_->residual(u, lambda, r, g);
    return std::hypot(p_.residual_norm(r), g.norm());
}
RunStatistics TimeIntegrator::advance(Vec &solution, bool forbid_allocations) {
    require(solution.size() == p_.layout.size() && p_.admissible(solution),
            "Invalid initial solution");
    const auto &c = p_.config;
    require(c.steps > 0 || c.final_time < p_.example.shock_time(),
            "Final time reaches simple-wave shock");
    RunStatistics stats;
    State initial_integral = p_.integral(solution);
    auto begin = std::chrono::steady_clock::now();
    AllocationGuard guard(forbid_allocations);
    while (c.steps > 0 ? stats.steps < c.steps : stats.final_time < c.final_time) {
        double rate = p_.freeze(solution);
        // Preserve the reference project's CFL definition: CFL=dt*rate.
        double step = c.cfl / rate;
        if (c.steps == 0) {
            double remaining = c.final_time - stats.final_time;
            step = std::min(step, remaining);
            if (remaining - step <=
                16 * std::numeric_limits<double>::epsilon() * std::max(remaining, c.final_time))
                step = remaining;
        }
        require(std::isfinite(step) && step > 0 && stats.final_time + step > stats.final_time,
                "Invalid CFL step");
        require(stats.final_time + step < p_.example.shock_time(),
                "Step reaches simple-wave shock");
        p_.step_size = step;
        state_ = solution;
        p_.trace(state_, lambda_);
        bool converged = false;
        for (int it = 0; it < c.max_newton; ++it) {
            double merit = evaluate(state_, lambda_, residual_, skeleton_);
            if (merit <= c.nonlinear_tolerance) {
                converged = true;
                break;
            }
            LinearResult linear{};
            if (method_ == Method::DG) {
                dg_->refresh(state_);
                rhs_ = -residual_;
                auto blocks = [&](int row, int col, BlockRef b) {
                    dg_->block(row, col, b);
                };
                if (direct_) {
                    direct_->refresh(blocks);
                    linear = direct_->solve(rhs_, direction_, linear_weights_);
                } else {
                    preconditioner_->refresh(blocks);
                    linear = gmres_->solve(Action::bind(*dg_), Action::bind(*preconditioner_), rhs_,
                                           direction_, c.max_linear, c.linear_tolerance);
                }
            } else {
                hybrid_->refresh(state_, lambda_);
                rhs_ = hybrid_->rhs();
                auto blocks = [&](int row, int col, BlockRef b) { hybrid_->block(row, col, b); };
                if (direct_) {
                    direct_->refresh(blocks);
                    linear = direct_->solve(rhs_, trace_direction_, linear_weights_);
                } else {
                    preconditioner_->refresh(blocks);
                    linear =
                        gmres_->solve(Action::bind(*hybrid_), Action::bind(*preconditioner_), rhs_,
                                      trace_direction_, c.max_linear, c.linear_tolerance);
                }
                hybrid_->recover(trace_direction_, direction_);
            }
            stats.max_linear_residual =
                std::max(stats.max_linear_residual, linear.relative_residual);
            stats.krylov_iterations += linear.iterations;
            ++stats.newton_iterations;
            bool accepted = false;
            double alpha = 1;
            for (int search = 0; search < 30; ++search) {
                trial_ = state_ + alpha * direction_;
                trial_lambda_ = lambda_;
                if (method_ == Method::HDG)
                    trial_lambda_.noalias() += alpha * trace_direction_;
                if (p_.admissible(trial_)) {
                    double next = evaluate(trial_, trial_lambda_, trial_residual_, trial_skeleton_);
                    if (next <= (1 - 1e-4 * alpha) * merit) {
                        accepted = true;
                        break;
                    }
                }
                alpha *= .5;
            }
            if (!accepted)
                throw std::runtime_error(
                    "Newton line search failed at time=" + std::to_string(stats.final_time) +
                    ", iteration=" + std::to_string(it) + ", residual=" + std::to_string(merit));
            state_.swap(trial_);
            lambda_.swap(trial_lambda_);
            if (evaluate(state_, lambda_, residual_, skeleton_) <= c.nonlinear_tolerance) {
                converged = true;
                break;
            }
        }
        require(converged, "Newton iteration limit reached");
        solution = state_;
        stats.final_time = c.steps == 0 && step >= c.final_time - stats.final_time
                               ? c.final_time
                               : stats.final_time + step;
        stats.minimum_step = stats.steps == 0 ? step : std::min(stats.minimum_step, step);
        stats.maximum_step = std::max(stats.maximum_step, step);
        ++stats.steps;
        stats.conservation_drift =
            std::max(stats.conservation_drift,
                     (p_.integral(solution) - initial_integral).cwiseAbs().maxCoeff());
    }
    stats.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    return stats;
}
} // namespace euler

#include "euler/time_integrator.hpp"
#include <algorithm>
#include <chrono>
#include <limits>
namespace euler {
namespace {
class AllocationGuard {
    bool previous_;
  public:
    explicit AllocationGuard(bool forbid) : previous_(Eigen::internal::is_malloc_allowed()) {
        if (forbid) Eigen::internal::set_is_malloc_allowed(false);
    }
    ~AllocationGuard() { Eigen::internal::set_is_malloc_allowed(previous_); }
};
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
} // namespace
TimeIntegrator::TimeIntegrator(Problem &p, Method method)
    : p_(p), method_(method), state_(p.layout.size()), trial_(p.layout.size()),
      residual_(p.layout.size()), trial_residual_(p.layout.size()), direction_(p.layout.size()),
      rhs_(method == Method::DG ? p.layout.size() : p.trace_size()),
      lambda_(p.trace_size()), trace_direction_(p.trace_size()),
      equation_error_(p.layout.size()), independent_error_(p.layout.size()),
      step_state_(p.layout.size()), stage_base_(p.layout.size()), slopes_(p.layout.size(), 4),
      linear_weights_(method == Method::DG ? p.inverse_mass : p.inverse_trace_mass) {
    if (method == Method::HDG) hybrid_.emplace(p);
    else dg_.emplace(p);
    if (method == Method::HDG && p.config.check_update) diagnostic_.emplace(p);
    const int width=method==Method::DG?p.layout.block_size():p.face_width();
    const int blocks=method==Method::DG?p.layout.cells:p.face_count();
    if(p.config.linear==LinearSolver::Direct) {
        require(rhs_.size()<=5000,"Direct reference limited to 5000 unknowns");
        const auto graph=method==Method::DG?BlockGraph::cartesian(p.layout.nx,p.layout.ny):p.trace_graph();
        direct_.emplace(graph,width);
    } else {
        preconditioner_.emplace(blocks,width);
        gmres_.emplace(rhs_.size(),p.config.restart,linear_weights_);
    }
}
double TimeIntegrator::dual_norm(const Vec &x) const {
    return std::sqrt((x.array().square() * p_.inverse_mass.array()).sum());
}
double TimeIntegrator::state_norm(const Vec &x) const {
    return std::sqrt((x.array().square() * p_.mass.array()).sum());
}
RunStatistics TimeIntegrator::advance(Vec &solution, bool forbid_allocations) {
    require(solution.size() == p_.layout.size() && p_.admissible(solution), "Invalid initial solution");
    const auto &c = p_.config;
    require(c.steps > 0 || c.final_time < p_.example.shock_time(), "Final time reaches simple-wave shock");
    constexpr double alpha = .24169426078821, beta = alpha / 4., eta = .12915286960590;
    constexpr double a[4][4] = {{alpha,0,0,0}, {-alpha,alpha,0,0}, {0,1-alpha,alpha,0},
                               {beta,eta,.5-beta-eta-alpha,alpha}};
    constexpr double b[4] = {0,1./6,1./6,2./3}, stage_c[4] = {alpha,0,1,.5};
    const int stages = c.integrator == "be" ? 1 : 4;
    RunStatistics stats;
    const double initial_rate = p_.freeze(solution);
    // Reserve outside the measured interval. Only standard-library storage can
    // grow if an unusually long run exceeds this estimate; Eigen stays fixed.
    const std::size_t estimate = c.steps > 0 ? c.steps :
        static_cast<std::size_t>(std::ceil(c.final_time * initial_rate / c.cfl) * 4 + 8);
    stats.updates.reserve(estimate * stages * (c.newton_steps > 0 ? c.newton_steps : c.max_newton));
    const State initial_integral = p_.integral(solution);
    auto begin = Clock::now();
    AllocationGuard guard(forbid_allocations);
    while (c.steps > 0 ? stats.steps < c.steps : stats.final_time < c.final_time) {
        const double rate = p_.freeze(solution);
        double step = c.cfl / rate;
        if (c.steps == 0) {
            const double remaining = c.final_time - stats.final_time;
            step = std::min(step, remaining);
            if (remaining - step <= 16 * std::numeric_limits<double>::epsilon() *
                                        std::max(remaining, c.final_time)) step = remaining;
        }
        require(std::isfinite(step) && step > 0 && stats.final_time + step > stats.final_time,
                "Invalid CFL step");
        require(stats.final_time + step < p_.example.shock_time(), "Step reaches simple-wave shock");
        step_state_ = solution;
        state_ = solution;
        for (int stage = 0; stage < stages; ++stage) {
            const double diagonal = stages == 1 ? 1. : a[stage][stage];
            stage_base_ = step_state_;
            for (int j = 0; j < stage; ++j) stage_base_.noalias() += step * a[stage][j] * slopes_.col(j);
            if (stage > 0) p_.freeze(state_);
            p_.step_size = diagonal * step;
            p_.old = stage_base_;
            const double tau = p_.step_size, old_norm = state_norm(stage_base_);
            const int limit = c.newton_steps > 0 ? c.newton_steps : c.max_newton;
            double initial_rms = -1.;
            bool converged = false, residual_ready = false;
            for (int it = 0; it <= limit; ++it) {
                const auto update_start = Clock::now();
                if (!residual_ready) p_.dg_residual(state_, residual_);
                const double norm = dual_norm(residual_), rms = norm / std::sqrt(p_.example.lx*p_.example.ly);
                if (initial_rms < 0) initial_rms = rms;
                if (c.newton_steps == 0 && rms <= c.nonlinear_tolerance * std::max(1., initial_rms)) {
                    converged = true;
                    break;
                }
                if (it == limit) break;
                // Convert the reference's mass-inverted stage equation allowance
                // to our physical weak equation, which is divided by tau.
                const double allowance = 128 * std::numeric_limits<double>::epsilon() *
                    ((state_norm(state_) + old_norm) / tau + norm);
                const double target = std::max(10 * c.linear_tolerance * norm, allowance);
                LinearResult linear{};
                if (method_ == Method::DG) {
                    dg_->refresh(state_);
                    rhs_ = -residual_;
                    auto blocks = [&](int row, int col, BlockRef block) { dg_->block(row, col, block); };
                    if (direct_) {
                        direct_->refresh(blocks);
                        linear = direct_->solve(rhs_, direction_, linear_weights_);
                        dg_->apply(direction_, equation_error_);
                        equation_error_ += residual_;
                    } else {
                        preconditioner_->refresh(blocks);
                        linear = gmres_->solve(Action::bind(*dg_), Action::bind(*preconditioner_),
                                               rhs_, direction_, c.max_linear, c.linear_tolerance);
                        equation_error_ = -gmres_->residual();
                    }
                } else {
                    hybrid_->refresh_projected(state_, residual_);
                    rhs_ = hybrid_->rhs();
                    auto blocks = [&](int row, int col, BlockRef block) { hybrid_->block(row, col, block); };
                    if (direct_) direct_->refresh(blocks);
                    else preconditioner_->refresh(blocks);
                    double tolerance = c.linear_tolerance;
                    bool valid = false;
                    for (int attempt = 0; attempt < 4; ++attempt) {
                        LinearResult current;
                        if (direct_) current = direct_->solve(rhs_, trace_direction_, linear_weights_);
                        else current = gmres_->solve(Action::bind(*hybrid_), Action::bind(*preconditioner_),
                                                      rhs_, trace_direction_, c.max_linear, tolerance);
                        linear.iterations += current.iterations;
                        linear.relative_residual = current.relative_residual;
                        hybrid_->recover(trace_direction_, direction_);
                        if(direct_) {
                            hybrid_->apply(trace_direction_,lambda_);
                            lambda_-=rhs_;
                        } else lambda_=-gmres_->residual();
                        hybrid_->recovered_residual_from_schur(lambda_,equation_error_);
                        valid = dual_norm(equation_error_) <= target;
                        if (valid) break;
                        require(!direct_, "Direct Schur solve fails DG residual criterion");
                        tolerance *= .01;
                    }
                    require(valid, "Recovered DG residual criterion failed after Schur refinement");
                }
                const double error_norm = dual_norm(equation_error_);
                require(error_norm <= target, "DG Newton linear residual too large");
                if (diagnostic_) {
                    diagnostic_->refresh(state_);
                    diagnostic_->apply(direction_, independent_error_);
                    independent_error_ += residual_;
                    require(dual_norm(independent_error_) <= std::max(4 * target, 1e-12),
                            "Independent recovered DG equation check failed");
                    independent_error_ -= equation_error_;
                    require(dual_norm(independent_error_) <= std::max(4 * allowance, 1e-12),
                            "Recovered residual identity failed");
                }
                double damping = 1.;
                bool accepted = false;
                for (int search = 0; search < 30; ++search) {
                    trial_ = state_ + damping * direction_;
                    if (p_.admissible(trial_)) {
                        p_.dg_residual(trial_, trial_residual_);
                        if (dual_norm(trial_residual_) <= (1 - 1e-4 * damping) * norm + allowance) {
                            accepted = true;
                            break;
                        }
                    }
                    damping *= .5;
                }
                require(accepted, "Newton line search failed");
                state_.swap(trial_);
                residual_.swap(trial_residual_);
                residual_ready = true;
                const double duration = elapsed(update_start);
                stats.max_linear_residual = std::max(stats.max_linear_residual, linear.relative_residual);
                const double dgrel = norm > 0 ? error_norm / norm : 0.;
                stats.max_dg_linear_residual = std::max(stats.max_dg_linear_residual, dgrel);
                stats.krylov_iterations += linear.iterations;
                if (stats.newton_iterations >= c.warmup_updates) {
                    stats.measured_newton_seconds += duration;
                    ++stats.measured_updates;
                }
                ++stats.newton_iterations;
                stats.updates.push_back({stats.steps, stage, it, linear.iterations,
                    stats.final_time + step, stats.final_time + (stages == 1 ? 1. : stage_c[stage]) * step,
                    step, step * rate, rms, damping, linear.relative_residual, dgrel, duration});
            }
            require(converged || c.newton_steps > 0, "Newton iteration limit reached");
            slopes_.col(stage) = (state_ - stage_base_) / tau;
        }
        if (stages == 1) solution = state_;
        else {
            solution = step_state_;
            for (int stage = 0; stage < stages; ++stage) solution.noalias() += step * b[stage] * slopes_.col(stage);
            require(p_.admissible(solution), "DIRK3 step produced an inadmissible state");
        }
        stats.final_time = c.steps == 0 && step >= c.final_time - stats.final_time
                               ? c.final_time : stats.final_time + step;
        stats.minimum_step = stats.steps == 0 ? step : std::min(stats.minimum_step, step);
        stats.maximum_step = std::max(stats.maximum_step, step);
        ++stats.steps;
        stats.conservation_drift = std::max(stats.conservation_drift,
            (p_.integral(solution) - initial_integral).cwiseAbs().maxCoeff());
    }
    stats.seconds = elapsed(begin);
    return stats;
}
} // namespace euler

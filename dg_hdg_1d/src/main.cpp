#include "euler/time_integrator.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
int main(int argc, char **argv) {
    using namespace euler;
    try {
        const Config config = parse(argc, argv);
        Problem problem(config);
        const auto prefix = std::filesystem::path(config.prefix);
        if (prefix.has_parent_path()) std::filesystem::create_directories(prefix.parent_path());
        std::ofstream summary(config.prefix + "_summary.csv"), updates(config.prefix + "_iterations.csv");
        require(bool(summary) && bool(updates), "Cannot open output files");
        summary << std::setprecision(17)
            << "method,cells,k,quad,cfl,integrator,linear,preconditioner,ordering,ntol,ltol,safety,"
               "fixed_newton_steps,warmup_updates,steps,final_time,newton_iterations,krylov_iterations,"
               "max_linear_relative_residual,max_dg_linear_relative_residual,min_step,max_step,mean_dt,"
               "time_loop_s,measured_updates,mean_newton_s,volume_dofs,trace_dofs,"
               "rho_L2,momentum_L2,energy_L2,max_conservation_drift\n";
        updates << std::setprecision(17)
            << "method,step,stage,iteration,time,stage_time,dt,cfl,residual_physical_rms,damping,"
               "krylov_iterations,linear_relative_residual,dg_linear_relative_residual,iteration_s,measured\n";
        std::ofstream solution_file;
        if (!config.timing_only) {
            solution_file.open(config.prefix + "_solution.csv");
            require(bool(solution_file), "Cannot open solution output");
            solution_file << std::setprecision(17)
                << "method,cell,node,x,rho,momentum,energy,rho_exact,momentum_exact,energy_exact\n";
        }
        Vec dg_solution(problem.layout.size()), solution(problem.layout.size());
        bool have_dg = false;
        for (Method method : {Method::DG, Method::HDG}) {
            const std::string name = method == Method::DG ? "dg" : "hdg";
            if (config.method != "both" && config.method != name) continue;
            problem.initial(solution);
            TimeIntegrator integrator(problem, method);
            const auto stats = integrator.advance(solution);
            const State error = config.timing_only ? State::Constant(std::numeric_limits<double>::quiet_NaN())
                                                  : problem.error(solution, stats.final_time);
            summary << name << ',' << config.cells << ',' << config.degree << ',' << problem.basis.points
                << ',' << config.cfl << ',' << config.integrator << ','
                << (config.linear == LinearSolver::Direct ? "direct" : "gmres") << ','
                << (config.linear == LinearSolver::Direct ? "none" :
                    config.preconditioner == Preconditioner::BlockJacobi ? "block-jacobi" : "block-ilu0") << ','
                << (config.linear == LinearSolver::Direct ? "none" :
                    config.preconditioner == Preconditioner::BlockJacobi || !config.colamd ? "natural" : "colamd")
                << ',' << config.nonlinear_tolerance << ',' << config.linear_tolerance << ',' << config.safety
                << ',' << config.newton_steps << ',' << config.warmup_updates << ',' << stats.steps
                << ',' << stats.final_time << ',' << stats.newton_iterations << ',' << stats.krylov_iterations
                << ',' << stats.max_linear_residual << ',' << stats.max_dg_linear_residual
                << ',' << stats.minimum_step << ',' << stats.maximum_step << ','
                << (stats.steps ? stats.final_time / stats.steps : 0.) << ',' << stats.seconds
                << ',' << stats.measured_updates << ','
                << (stats.measured_updates ? stats.measured_newton_seconds / stats.measured_updates : 0.)
                << ',' << problem.layout.size() << ',' << problem.layout.cells * NVAR;
            for (int v = 0; v < NVAR; ++v) summary << ',' << error[v];
            summary << ',' << stats.conservation_drift << '\n';
            for (std::size_t i = 0; i < stats.updates.size(); ++i) {
                const auto &u = stats.updates[i];
                updates << name << ',' << u.step << ',' << u.stage << ',' << u.iteration << ','
                    << u.time << ',' << u.stage_time << ',' << u.dt << ',' << u.cfl << ','
                    << u.residual_physical_rms << ',' << u.damping << ',' << u.krylov_iterations << ','
                    << u.linear_relative_residual << ',' << u.dg_linear_relative_residual << ','
                    << u.iteration_s << ',' << (i >= static_cast<std::size_t>(config.warmup_updates)) << '\n';
            }
            if (!config.timing_only)
                for (int cell = 0; cell < problem.layout.cells; ++cell)
                    for (int i = 0; i < problem.layout.nodes; ++i) {
                        const double x = -1 + problem.h * (cell + .5) + .5 * problem.h * problem.basis.x[i];
                        const State exact = problem.example.exact(x, stats.final_time);
                        solution_file << name << ',' << cell << ',' << i << ',' << x;
                        for (int v = 0; v < NVAR; ++v) solution_file << ',' << solution[problem.layout.index(cell,v,i)];
                        for (int v = 0; v < NVAR; ++v) solution_file << ',' << exact[v];
                        solution_file << '\n';
                    }
            std::cout << name << ": steps=" << stats.steps << ", Newton=" << stats.newton_iterations
                << ", Krylov=" << stats.krylov_iterations << ", time=" << stats.seconds
                << " s, mean Newton=" << (stats.measured_updates ? stats.measured_newton_seconds / stats.measured_updates : 0.)
                << " s, conservation drift=" << stats.conservation_drift << '\n';
            if (method == Method::DG) { dg_solution = solution; have_dg = true; }
            else if (have_dg) {
                const Vec difference = dg_solution - solution;
                const double maximum = difference.cwiseAbs().maxCoeff();
                std::ofstream comparison(config.prefix + "_difference.csv");
                require(bool(comparison), "Cannot open comparison output");
                comparison << std::setprecision(17) << "rho_L2,momentum_L2,energy_L2,joint_L2,max_nodal\n";
                double joint = 0.;
                for (int v = 0; v < NVAR; ++v) {
                    double squared = 0.;
                    for (int cell = 0; cell < problem.layout.cells; ++cell)
                        squared += (problem.layout.cell(difference,cell).col(v).array().square() *
                                    problem.layout.cell(problem.mass,cell).col(v).array()).sum();
                    comparison << std::sqrt(squared) << ',';
                    joint += squared;
                }
                comparison << std::sqrt(joint) << ',' << maximum << '\n';
                std::cout << "DG/HDG max difference=" << maximum << '\n';
                require(maximum < 1e-8, "DG/HDG solutions disagree");
            }
        }
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

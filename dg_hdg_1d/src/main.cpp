#include "euler/time_integrator.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
    using namespace euler;
    try {
        Config config = parse(argc, argv);
        Problem problem(config);
        auto prefix = std::filesystem::path(config.prefix);
        if (prefix.has_parent_path())
            std::filesystem::create_directories(prefix.parent_path());
        std::ofstream summary(config.prefix + "_summary.csv"),
            solution_file(config.prefix + "_solution.csv");
        require(bool(summary) && bool(solution_file), "Cannot open output files");
        summary << std::setprecision(17)
                << "method,cells,k,quad,cfl,linear,preconditioner,ordering,steps,final_time,newton_"
                   "iterations,krylov_iterations,max_linear_relative_residual,min_step,max_step,"
                   "time_loop_s,rho_L2,momentum_L2,energy_L2,max_conservation_drift\n";
        solution_file
            << std::setprecision(17)
            << "method,cell,node,x,rho,momentum,energy,rho_exact,momentum_exact,energy_exact\n";
        Vec dg_solution(problem.layout.size()), solution(problem.layout.size());
        bool have_dg = false;
        for (Method method : {Method::DG, Method::HDG}) {
            std::string name = method == Method::DG ? "dg" : "hdg";
            if (config.method != "both" && config.method != name)
                continue;
            problem.initial(solution);
            TimeIntegrator integrator(problem, method);
            auto stats = integrator.advance(solution);
            State error = problem.error(solution, stats.final_time);
            summary << name << ',' << config.cells << ',' << config.degree << ','
                    << problem.basis.points << ',' << config.cfl << ','
                    << (config.linear == LinearSolver::Direct ? "direct" : "gmres") << ','
                    << (config.linear == LinearSolver::Direct                ? "none"
                        : config.preconditioner == Preconditioner::BlockILU0 ? "block-ilu0"
                                                                             : "block-jacobi")
                    << ','
                    << (config.linear == LinearSolver::Direct ? "none"
                        : config.colamd                       ? "colamd"
                                                              : "natural")
                    << ',' << stats.steps << ',' << stats.final_time << ','
                    << stats.newton_iterations << ',' << stats.krylov_iterations << ','
                    << stats.max_linear_residual << ',' << stats.minimum_step << ','
                    << stats.maximum_step << ',' << stats.seconds;
            for (int v = 0; v < NVAR; ++v)
                summary << ',' << error[v];
            summary << ',' << stats.conservation_drift << '\n';
            for (int cell = 0; cell < problem.layout.cells; ++cell)
                for (int i = 0; i < problem.layout.nodes; ++i) {
                    double x = -1 + problem.h * (cell + .5) + .5 * problem.h * problem.basis.x[i];
                    State exact = problem.example.exact(x, stats.final_time);
                    solution_file << name << ',' << cell << ',' << i << ',' << x;
                    for (int v = 0; v < NVAR; ++v)
                        solution_file << ',' << solution[problem.layout.index(cell, v, i)];
                    for (int v = 0; v < NVAR; ++v)
                        solution_file << ',' << exact[v];
                    solution_file << '\n';
                }
            std::cout << name << ": steps=" << stats.steps << ", Newton=" << stats.newton_iterations
                      << ", Krylov=" << stats.krylov_iterations << ", time=" << stats.seconds
                      << " s, conservation drift=" << stats.conservation_drift << '\n';
            if (method == Method::DG) {
                dg_solution = solution;
                have_dg = true;
            } else if (have_dg) {
                double difference = (dg_solution - solution).cwiseAbs().maxCoeff();
                std::cout << "DG/HDG max difference=" << difference << '\n';
                require(difference < 1e-8, "DG/HDG solutions disagree");
            }
        }
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

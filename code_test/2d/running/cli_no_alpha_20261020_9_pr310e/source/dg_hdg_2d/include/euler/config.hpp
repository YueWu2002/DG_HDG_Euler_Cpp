#pragma once
#include <string>
namespace euler {
enum class Method { DG, HDG };
enum class Preconditioner { BlockJacobi, BlockILU0 };
enum class LinearSolver { GMRES, Direct };
struct Config {
    int nx = 8, ny = 8, degree = 3, quadrature = 0;
    int restart = 20, max_linear = 1500, max_newton = 30, steps = 0;
    int newton_steps = 0, warmup_updates = 0;
    double cfl = 2., final_time = .2, safety = 1.1;
    double nonlinear_tolerance = 1e-12, linear_tolerance = 1e-12;
    Preconditioner preconditioner = Preconditioner::BlockJacobi;
    LinearSolver linear = LinearSolver::GMRES;
    bool colamd = false, timing_only = false, check_update = false, save_solution = true;
    std::string integrator = "ssp-dirk3";
    std::string method = "both", testcase = "diagonal", prefix = "output/run_2d";
};
Config parse(int argc, char **argv);
void validate(const Config &config);
} // namespace euler

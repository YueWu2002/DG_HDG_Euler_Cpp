#pragma once
#include <string>
namespace euler {
enum class Method { DG, HDG };
enum class Preconditioner { BlockJacobi, BlockILU0 };
enum class LinearSolver { GMRES, Direct };
struct Config {
    int nx = 8, ny = 8, degree = 3, quadrature = 0;
    int restart = 30, max_linear = 600, max_newton = 25, steps = 0;
    double cfl = .2, final_time = .1, safety = 1.1;
    double nonlinear_tolerance = 1e-11, linear_tolerance = 1e-12;
    Preconditioner preconditioner = Preconditioner::BlockILU0;
    LinearSolver linear = LinearSolver::GMRES;
    bool colamd = true;
    std::string method = "both", testcase = "diagonal", prefix = "output/run_2d";
};
Config parse(int argc, char **argv);
void validate(const Config &config);
} // namespace euler

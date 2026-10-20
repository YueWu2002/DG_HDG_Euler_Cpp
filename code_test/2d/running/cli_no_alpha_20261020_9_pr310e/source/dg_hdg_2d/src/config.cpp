#include "euler/config.hpp"
#include "euler/types.hpp"
#include <cstdlib>
#include <iostream>
namespace euler {
void validate(const Config &c) {
    require(c.nx >= 3 && c.ny >= 3 && c.degree >= 0 && c.degree <= 8 && c.quadrature >= 0, "Invalid mesh or degree");
    require(c.testcase == "diagonal" || c.testcase == "wave" || c.testcase == "axis", "Invalid testcase");
    require(c.quadrature == 0 || c.quadrature >= c.degree + 1, "Quadrature must be >= degree+1");
    require(c.restart > 0 && c.max_linear > 0 && c.max_newton > 0 && c.steps >= 0,
            "Invalid iteration limit");
    require(c.newton_steps >= 0 && c.warmup_updates >= 0, "Invalid timing update count");
    require(c.integrator == "be" || c.integrator == "ssp-dirk3", "Invalid integrator");
    require(std::isfinite(c.cfl) && c.cfl > 0 && std::isfinite(c.final_time) && c.final_time >= 0,
            "CFL must be positive and final time nonnegative");
    require(std::isfinite(c.safety) && c.safety >= 1 && std::isfinite(c.nonlinear_tolerance) &&
                c.nonlinear_tolerance > 0 && std::isfinite(c.linear_tolerance) &&
                c.linear_tolerance > 0 && c.linear_tolerance < 1,
            "Invalid solver tolerance or wave-speed safety");
    require(c.method == "dg" || c.method == "hdg" || c.method == "both",
            "Method must be dg, hdg, or both");
}
Config parse(int argc, char **argv) {
    Config c;
    bool time_seen = false, steps_seen = false;
    for (int i = 1; i < argc; ++i) {
        std::string option = argv[i];
        auto value = [&]() {
            require(i + 1 < argc, "Missing option value");
            return std::string(argv[++i]);
        };
        auto integer = [&]() {
            auto text = value();
            std::size_t end;
            int n = std::stoi(text, &end);
            require(end == text.size(), "Invalid integer");
            return n;
        };
        auto real = [&]() {
            auto text = value();
            std::size_t end;
            double n = std::stod(text, &end);
            require(end == text.size(), "Invalid number");
            return n;
        };
        if (option == "--nx")
            c.nx = integer();
        else if (option == "--ny")
            c.ny = integer();
        else if (option == "--testcase")
            c.testcase = value();
        else if (option == "--no-save-solution")
            c.save_solution = false;
        else if (option == "--newton-steps")
            c.newton_steps = integer();
        else if (option == "--warmup-updates")
            c.warmup_updates = integer();
        else if (option == "--integrator")
            c.integrator = value();
        else if (option == "--timing-only")
            c.timing_only = true;
        else if (option == "--check-update")
            c.check_update = true;
        else if (option == "--k")
            c.degree = integer();
        else if (option == "--quad")
            c.quadrature = integer();
        else if (option == "--cfl")
            c.cfl = real();
        else if (option == "--T") {
            c.final_time = real();
            time_seen = true;
        } else if (option == "--steps") {
            c.steps = integer();
            steps_seen = true;
        } else if (option == "--safety")
            c.safety = real();
        else if (option == "--method")
            c.method = value();
        else if (option == "--prefix")
            c.prefix = value();
        else if (option == "--ntol")
            c.nonlinear_tolerance = real();
        else if (option == "--ltol")
            c.linear_tolerance = real();
        else if (option == "--restart")
            c.restart = integer();
        else if (option == "--max-linear")
            c.max_linear = integer();
        else if (option == "--max-newton")
            c.max_newton = integer();
        else if (option == "--linear") {
            auto name = value();
            require(name == "direct" || name == "gmres", "Invalid linear solver");
            c.linear = name == "direct" ? LinearSolver::Direct : LinearSolver::GMRES;
        } else if (option == "--ordering") {
            auto order = value();
            require(order == "natural" || order == "colamd", "Invalid ordering");
            c.colamd = order == "colamd";
        } else if (option == "--precond") {
            auto name = value();
            require(name == "block-jacobi", "2D currently supports block-jacobi only");
            c.preconditioner =
                name == "block-jacobi" ? Preconditioner::BlockJacobi : Preconditioner::BlockILU0;
        } else if (option == "--help") {
            std::cout
                << "2D periodic tensor-product Euler DG / condensed HDG\n"
                   "--nx NX --ny NY --testcase diagonal|wave|axis --k K --quad Q --cfl C (--T T | --steps N)\n"
                   "--method dg|hdg|both --linear direct|gmres --precond block-jacobi\n"
                   "--ordering natural|colamd --restart N --max-linear N --max-newton N\n"
                   "--ntol R --ltol R --safety S --prefix PATH\n"
                   "--integrator ssp-dirk3|be --newton-steps N --warmup-updates N\n"
                   "--timing-only --check-update --no-save-solution\n";
            std::exit(0);
        } else
            throw std::runtime_error("Unknown option: " + option);
    }
    require(!(time_seen && steps_seen), "Choose either --T or --steps");
    require(!steps_seen || c.steps > 0, "steps must be positive");
    validate(c);
    return c;
}
} // namespace euler

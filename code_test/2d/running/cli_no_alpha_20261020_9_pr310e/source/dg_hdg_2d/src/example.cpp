#include "euler/example.hpp"
#include <limits>
namespace euler {
double Example::shock_time() const {
    if (testcase == "wave")
        return std::numeric_limits<double>::infinity();
    return 4 / ((gas.gamma() + 1) * std::acos(-1.));
}
State Example::exact(double x, double y, double t) const {
    double g = gas.gamma(), pi = std::acos(-1.);
    State u;
    if (testcase == "wave") {
        double rho = 1 + .2 * std::sin(pi * (x - t)) * std::sin(pi * (y - .5 * t));
        u << rho, rho, .5 * rho, 1 / (g - 1) + .625 * rho;
        return u;
    }
    require(t < shock_time(), "Simple wave reaches shock time");
    bool diagonal = testcase == "diagonal";
    if (diagonal) {
        x = cosine * x + sine * y;
    }
    double a = (3 - g) / 2, b = (g + 1) / 4;
    double lo = x - t * (a + b), hi = x - t * (a - b), xi = x - t * a;
    for (int it = 0; it < 80; ++it) {
        double f = xi + t * (a + b * std::sin(pi * xi)) - x;
        if (std::abs(f) < 2e-15)
            break;
        if (f > 0) hi = xi; else lo = xi;
        double z = xi - f / (1 + t * b * pi * std::cos(pi * xi));
        xi = z > lo && z < hi ? z : .5 * (lo + hi);
    }
    double jm = std::sin(pi * xi), velocity = .5 * (2 + jm), sound = (g - 1) * (2 - jm) / 4;
    double rho = std::pow(sound * sound / g, 1 / (g - 1));
    u[Density] = rho;
    u[MomentumX] = rho * velocity * (diagonal ? cosine : 1.);
    u[MomentumY] = diagonal ? rho * velocity * sine : 0.;
    u[Energy] = std::pow(rho, g) / (g - 1) + .5 * rho * velocity * velocity;
    return u;
}
} // namespace euler

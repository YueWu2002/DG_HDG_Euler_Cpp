#include "euler/simple_wave.hpp"
namespace euler {
State SimpleWave::exact(double x, double t) const {
    require(t < shock_time(), "Simple wave reaches shock time");
    double g = gas.gamma(), pi = std::acos(-1.0), a = (3 - g) / 2, b = (g + 1) / 4;
    double lo = x - t * (a + b), hi = x - t * (a - b), xi = x - t * a;
    for (int i = 0; i < 80; ++i) {
        double f = xi + t * (a + b * std::sin(pi * xi)) - x;
        if (std::abs(f) < 2e-15)
            break;
        if (f > 0)
            hi = xi;
        else
            lo = xi;
        double z = xi - f / (1 + t * b * pi * std::cos(pi * xi));
        xi = z > lo && z < hi ? z : .5 * (lo + hi);
    }
    double jm = std::sin(pi * xi), v = .5 * (2 + jm), c = (g - 1) * (2 - jm) / 4;
    double rho = std::pow(c * c / g, 1 / (g - 1)), p = std::pow(rho, g);
    State u;
    u[Density] = rho;
    u[Momentum] = rho * v;
    u[Energy] = p / (g - 1) + .5 * rho * v * v;
    return u;
}
} // namespace euler

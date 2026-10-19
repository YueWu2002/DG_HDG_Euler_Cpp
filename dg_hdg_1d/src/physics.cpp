#include "euler/physics.hpp"
#include <algorithm>
namespace euler {
double Gas::pressure(const State &u) const {
    return (gamma_ - 1) * (u[Energy] - .5 * u[Momentum] * u[Momentum] / u[Density]);
}
State Gas::flux(const State &u) const {
    double velocity = u[Momentum] / u[Density], p = pressure(u);
    State f;
    f[Density] = u[Momentum];
    f[Momentum] = u[Momentum] * velocity + p;
    f[Energy] = (u[Energy] + p) * velocity;
    return f;
}
FluxJacobian Gas::jacobian(const State &u) const {
    double v = u[Momentum] / u[Density], H = (u[Energy] + pressure(u)) / u[Density], g = gamma_;
    FluxJacobian a;
    a << 0, 1, 0, .5 * (g - 3) * v * v, (3 - g) * v, g - 1, v * (.5 * (g - 1) * v * v - H),
        H - (g - 1) * v * v, g * v;
    return a;
}
double Gas::spectral_radius(const State &u) const {
    return std::abs(u[Momentum] / u[Density]) + std::sqrt(gamma_ * pressure(u) / u[Density]);
}
Speed einfeldt(const State &l, const State &r, const Gas &gas, double safety) {
    require(l[Density] > 0 && r[Density] > 0 && gas.pressure(l) > 0 && gas.pressure(r) > 0,
            "Invalid frozen face state");
    double sl = std::sqrt(l[Density]), sr = std::sqrt(r[Density]);
    double ul = l[Momentum] / l[Density], ur = r[Momentum] / r[Density];
    double cl = std::sqrt(gas.gamma() * gas.pressure(l) / l[Density]);
    double cr = std::sqrt(gas.gamma() * gas.pressure(r) / r[Density]);
    double v = (sl * ul + sr * ur) / (sl + sr);
    double H = (sl * (l[Energy] + gas.pressure(l)) / l[Density] +
                sr * (r[Energy] + gas.pressure(r)) / r[Density]) /
               (sl + sr);
    double c = std::sqrt((gas.gamma() - 1) * (H - .5 * v * v));
    return {safety * std::max({0.0, -ul + cl, -v + c}), safety * std::max({0.0, ur + cr, v + c})};
}
} // namespace euler

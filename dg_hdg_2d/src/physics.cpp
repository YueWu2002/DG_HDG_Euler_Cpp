#include "euler/physics.hpp"
#include <algorithm>
namespace euler {
double Gas::pressure(const State &u) const {
    return (gamma_ - 1) * (u[Energy] - .5 * u.segment<NDIM>(MomentumX).squaredNorm() / u[Density]);
}
Primitive Gas::primitive(const State &u) const {
    Primitive p;
    p.inverse_density = 1 / u[Density];
    p.velocity = u.segment<NDIM>(MomentumX) * p.inverse_density;
    p.pressure = (gamma_ - 1) * (u[Energy] - .5 * u.segment<NDIM>(MomentumX).dot(p.velocity));
    p.enthalpy = (u[Energy] + p.pressure) * p.inverse_density;
    return p;
}
State Gas::flux(const State &u, const Primitive &p, int axis) const {
    State f;
    f[Density] = u[MomentumX + axis];
    f.segment<NDIM>(MomentumX) = u.segment<NDIM>(MomentumX) * p.velocity[axis];
    f[MomentumX + axis] += p.pressure;
    f[Energy] = (u[Energy] + p.pressure) * p.velocity[axis];
    return f;
}
FluxJacobian Gas::jacobian(const State &u, const Primitive &p, int axis) const {
    State dp = State::Zero(), dv = State::Zero();
    dp[Density] = .5 * (gamma_ - 1) * p.velocity.squaredNorm();
    dp.segment<NDIM>(MomentumX) = -(gamma_ - 1) * p.velocity;
    dp[Energy] = gamma_ - 1;
    dv[Density] = -p.velocity[axis] * p.inverse_density;
    dv[MomentumX + axis] = p.inverse_density;
    FluxJacobian a = FluxJacobian::Zero();
    a(Density, MomentumX + axis) = 1;
    for (int d = 0; d < NDIM; ++d) {
        int momentum = MomentumX + d;
        a.row(momentum) = u[momentum] * dv.transpose();
        a(momentum, momentum) += p.velocity[axis];
    }
    a.row(MomentumX + axis) += dp.transpose();
    a.row(Energy) = (u[Energy] + p.pressure) * dv.transpose() + p.velocity[axis] * dp.transpose();
    a(Energy, Energy) += p.velocity[axis];
    return a;
}
Speed einfeldt(const State &l, const State &r, const Gas &gas, double safety, int axis) {
    auto pl = gas.primitive(l), pr = gas.primitive(r);
    require(l[Density] > 0 && r[Density] > 0 && pl.pressure > 0 && pr.pressure > 0,
            "Invalid frozen face state");
    double sl = std::sqrt(l[Density]), sr = std::sqrt(r[Density]);
    auto velocity = ((sl * pl.velocity + sr * pr.velocity) / (sl + sr)).eval();
    double H = (sl * pl.enthalpy + sr * pr.enthalpy) / (sl + sr);
    double sound = std::sqrt((gas.gamma() - 1) * (H - .5 * velocity.squaredNorm()));
    double cl = std::sqrt(gas.gamma() * pl.pressure * pl.inverse_density);
    double cr = std::sqrt(gas.gamma() * pr.pressure * pr.inverse_density);
    return {safety * std::max({0., -pl.velocity[axis] + cl, -velocity[axis] + sound}),
            safety * std::max({0., pr.velocity[axis] + cr, velocity[axis] + sound})};
}
} // namespace euler

#pragma once
#include "types.hpp"
namespace euler {
struct Primitive {
    double inverse_density, pressure, enthalpy;
    Eigen::Matrix<double, NDIM, 1> velocity;
};
class Gas {
  public:
    explicit Gas(double gamma) : gamma_(gamma) {}
    double gamma() const { return gamma_; }
    double pressure(const State &state) const;
    Primitive primitive(const State &state) const;
    State flux(const State &state, const Primitive &primitive, int axis) const;
    FluxJacobian jacobian(const State &state, const Primitive &primitive, int axis) const;
    State flux(const State &state, int axis) const { return flux(state, primitive(state), axis); }
    FluxJacobian jacobian(const State &state, int axis) const {
        return jacobian(state, primitive(state), axis);
    }

  private:
    double gamma_;
};
struct Speed { double a, b; };
Speed einfeldt(const State &left, const State &right, const Gas &gas, double safety, int axis);
} // namespace euler

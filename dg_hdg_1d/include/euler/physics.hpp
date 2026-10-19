#pragma once
#include "types.hpp"
namespace euler {
// Gas parameters belong to the equation/case, never to individual states.
class Gas {
  public:
    explicit Gas(double gamma) : gamma_(gamma) { require(gamma > 1, "gamma must exceed one"); }
    double gamma() const { return gamma_; }
    double pressure(const State &u) const;
    State flux(const State &u) const;
    FluxJacobian jacobian(const State &u) const;
    double spectral_radius(const State &u) const;

  private:
    double gamma_;
};
struct Speed {
    double a, b;
};
Speed einfeldt(const State &left, const State &right, const Gas &gas, double safety);
} // namespace euler

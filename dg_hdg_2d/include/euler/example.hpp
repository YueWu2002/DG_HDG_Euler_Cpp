#pragma once
#include "physics.hpp"
#include <string>
namespace euler {
struct Example {
    Gas gas{3.0};
    std::string testcase;
    double cosine, sine, lx, ly;
    explicit Example(std::string name, double degrees = 45.)
        : testcase(std::move(name)), cosine(std::cos(degrees * std::acos(-1.) / 180.)),
          sine(std::sin(degrees * std::acos(-1.) / 180.)),
          lx(testcase == "diagonal" ? 2. / cosine : 2.),
          ly(testcase == "diagonal" ? 2. / sine : 2.) {}
    State exact(double x, double y, double time) const;
    double shock_time() const;
};
} // namespace euler

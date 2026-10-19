#pragma once
#include "physics.hpp"
#include <string>
namespace euler {
struct Example {
    Gas gas{3.0};
    std::string testcase;
    explicit Example(std::string name) : testcase(std::move(name)) {}
    State exact(double x, double y, double time) const;
    double shock_time() const;
};
} // namespace euler

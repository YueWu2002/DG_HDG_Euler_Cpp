#pragma once
#include "euler/discretization.hpp"
namespace euler::check::reference {
// Independent general-quadrature implementation retained only for verification
// and small-grid kernel timing. It intentionally does not use cached coefficients.
void volume_block(const Problem &p, const Vec &u, int cell, BlockRef result);
void dg_action(const Problem &p, const Vec &u, const Vec &direction, Vec &result);
} // namespace euler::check::reference

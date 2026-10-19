#include "dg_reference.hpp"
#include "euler/dg_operator.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
namespace {
volatile double checksum = 0;
template <class Kernel> double microseconds(int repeats, Kernel &&kernel) {
    for (int i = 0; i < 10; ++i)
        kernel();
    auto begin = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i)
        kernel();
    return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - begin)
               .count() / repeats;
}
} // namespace
int main() {
    using namespace euler;
    std::cout << std::fixed << std::setprecision(3)
              << "cells,k,quad,reference_action_us,cached_action_us,refresh_us,"
                 "reference_volume_blocks_us,cached_volume_blocks_us,"
                 "reference_20_actions_us,refresh_plus_20_actions_us\n";
    for (int degree : {0, 1, 3, 6})
        for (int extra : {0, 2}) {
            Config c;
            c.cells = 32;
            c.degree = degree;
            c.quadrature = degree + 1 + extra;
            Problem p(c);
            Vec u(p.layout.size()), x = Vec::LinSpaced(p.layout.size(), -.2, .3), y(x.size());
            p.initial(u);
            p.freeze(u);
            p.step_size = .003;
            DGOperator op(p);
            VolumeJacobian volume(p);
            Mat block(p.layout.block_size(), p.layout.block_size());
            op.refresh(u);
            volume.refresh(u);
            Eigen::internal::set_is_malloc_allowed(false);
            auto reference_action = [&]() {
                check::reference::dg_action(p, u, x, y);
                checksum = y[0];
            };
            auto cached_action = [&]() {
                op.apply(x, y);
                checksum = y[0];
            };
            double old_action = microseconds(1000, reference_action);
            double new_action = microseconds(1000, cached_action);
            double refresh = microseconds(1000, [&]() { op.refresh(u); });
            double old_blocks = microseconds(300, [&]() {
                for (int cell = 0; cell < c.cells; ++cell) {
                    check::reference::volume_block(p, u, cell, block);
                    checksum = block(0, 0);
                }
            });
            double new_blocks = microseconds(300, [&]() {
                for (int cell = 0; cell < c.cells; ++cell) {
                    volume.block(cell, block);
                    checksum = block(0, 0);
                }
            });
            double old_cycle = microseconds(100, [&]() {
                for (int apply = 0; apply < 20; ++apply)
                    reference_action();
            });
            double new_cycle = microseconds(100, [&]() {
                op.refresh(u);
                for (int apply = 0; apply < 20; ++apply)
                    cached_action();
            });
            Eigen::internal::set_is_malloc_allowed(true);
            std::cout << c.cells << ',' << degree << ',' << p.basis.points << ',' << old_action
                      << ',' << new_action << ',' << refresh << ',' << old_blocks << ','
                      << new_blocks << ',' << old_cycle << ',' << new_cycle << '\n';
        }
}

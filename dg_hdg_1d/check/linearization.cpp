#include "check_1d.hpp"
#include "dg_reference.hpp"
#include "euler/dg_operator.hpp"
#include <iostream>
namespace euler::check {
namespace {
class NoEigenAllocations {
    bool previous_ = Eigen::internal::is_malloc_allowed();

  public:
    NoEigenAllocations() { Eigen::internal::set_is_malloc_allowed(false); }
    ~NoEigenAllocations() { Eigen::internal::set_is_malloc_allowed(previous_); }
};
} // namespace
void linearization_lifecycle() {
    for (int degree : {0, 1, 3, 6})
        for (int extra : {0, 2}) {
            Config c;
            c.cells = 5;
            c.degree = degree;
            c.quadrature = degree + 1 + extra;
            Problem p(c);
            Vec state(p.layout.size()), direction = Vec::LinSpaced(p.layout.size(), -.1, .2);
            Vec old_action(direction.size()), cached(direction.size()), expected(direction.size());
            p.initial(state);
            p.freeze(state);
            p.step_size = .003;
            DGOperator op(p);
            VolumeJacobian volume(p);
            const int width = p.layout.block_size();
            Mat old_block(width, width), block(width, width), expected_volume(width, width);
            {
                // No warm-up: check the very first coefficient refresh and Jv/block calls.
                NoEigenAllocations guard;
                op.refresh(state);
                volume.refresh(state);
                op.apply(direction, old_action);
                op.block(0, 0, old_block);
                for (int cell = 0; cell < c.cells; ++cell) {
                    volume.block(cell, block);
                    reference::volume_block(p, state, cell, expected_volume);
                    require((block - expected_volume).norm() <= 1e-13 * expected_volume.norm(),
                            "Collocated/overintegrated volume block reference mismatch");
                }
            }
            // Reusing the same state vector's memory must not silently change the snapshot.
            state += .0001 * direction;
            p.step_size *= .7;
            p.freeze(state);
            {
                NoEigenAllocations guard;
                op.apply(direction, cached);
                op.block(0, 0, block);
                require((cached - old_action).norm() == 0 && (block - old_block).norm() == 0,
                        "DG snapshot changed without refresh");
                op.refresh(state);
                op.apply(direction, cached);
                reference::dg_action(p, state, direction, expected);
                require((cached - expected).norm() <= 1e-13 * expected.norm(),
                        "DG refresh failed to update state, step, or frozen face speeds");
                require((cached - old_action).norm() > 1e-8 * old_action.norm(),
                        "DG refresh test did not change the linearization");
            }
        }
    std::cout << "allocation-free first refresh, coefficient reuse, snapshot invalidation passed\n";
}
} // namespace euler::check

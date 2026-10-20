#pragma once
#include "block_graph.hpp"
#include "dense_lu.hpp"
#include "euler/config.hpp"
namespace euler {
// ILU(0) means no fill in the permuted BLOCK graph. It is not scalar ILU.
// Matrix-free use: refresh asks for individual blocks; no scalar global matrix
// is assembled. An incomplete factorization necessarily stores block factors.
class BlockPreconditioner {
  public:
    BlockPreconditioner(const BlockGraph &graph, int width, Preconditioner kind, bool colamd,
                        int nx_nodes = 0, int ny_nodes = 1);
    template <class BlockProvider> void refresh(BlockProvider &&get_block) {
        for (int row = 0; row < factors_.graph.rows; ++row)
            for (int p = factors_.graph.row_ptr[row]; p < factors_.graph.row_ptr[row + 1]; ++p)
                get_block(order_[row], order_[factors_.graph.col_idx[p]], factors_.block(p));
        factorize();
    }
    void apply(const Vec &rhs, Vec &result) const;
    const BlockMatrix &factors() const { return factors_; }

  private:
    std::vector<int> order_;
    BlockMatrix factors_;
    std::vector<DenseLU> diagonal_;
    mutable Vec work_;
    Preconditioner kind_;
    // Symbolic update schedule, constructed once. No graph searches while
    // refactorizing or applying the preconditioner.
    std::vector<int> diagonal_slots_, update_offsets_, update_targets_, update_upper_;
    void factorize();
};
} // namespace euler

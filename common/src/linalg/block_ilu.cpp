#include "euler/linalg/block_ilu.hpp"
namespace euler {
static BlockGraph permuted_graph(const BlockGraph &graph, const std::vector<int> &order,
                                 Preconditioner kind) {
    std::vector<int> inverse(order.size());
    for (int i = 0; i < graph.rows; ++i)
        inverse[order[i]] = i;
    std::vector<std::vector<int>> rows(graph.rows);
    for (int i = 0; i < graph.rows; ++i)
        for (int p = graph.row_ptr[order[i]]; p < graph.row_ptr[order[i] + 1]; ++p)
            if (kind == Preconditioner::BlockILU0 || graph.col_idx[p] == order[i])
                rows[i].push_back(inverse[graph.col_idx[p]]);
    return BlockGraph(rows, graph.cells_x, graph.cells_y);
}
BlockPreconditioner::BlockPreconditioner(const BlockGraph &graph, int width, Preconditioner kind,
                                         bool colamd, int nx_nodes, int ny_nodes)
    : order_(block_ordering(graph, colamd)),
      factors_(permuted_graph(graph, order_, kind), width, nx_nodes, ny_nodes),
      work_(graph.rows * width), kind_(kind) {
    diagonal_.reserve(graph.rows);
    for (int row = 0; row < graph.rows; ++row)
        diagonal_.emplace_back(width);
    const auto &pattern = factors_.graph;
    diagonal_slots_.resize(pattern.rows);
    update_offsets_.resize(pattern.col_idx.size() + 1);
    for (int row = 0; row < pattern.rows; ++row) {
        diagonal_slots_[row] = pattern.find(row, row);
        for (int slot = pattern.row_ptr[row]; slot < pattern.row_ptr[row + 1]; ++slot) {
            update_offsets_[slot] = static_cast<int>(update_targets_.size());
            int eliminated_row = pattern.col_idx[slot];
            if (eliminated_row >= row)
                continue;
            for (int upper = pattern.row_ptr[eliminated_row];
                 upper < pattern.row_ptr[eliminated_row + 1]; ++upper) {
                int column = pattern.col_idx[upper];
                if (column <= eliminated_row)
                    continue;
                int target = pattern.find(row, column);
                if (target >= 0) {
                    update_targets_.push_back(target);
                    update_upper_.push_back(upper);
                }
            }
        }
    }
    update_offsets_.back() = static_cast<int>(update_targets_.size());
}
void BlockPreconditioner::factorize() {
    const auto &graph = factors_.graph;
    for (int i = 0; i < graph.rows; ++i) {
        if (kind_ == Preconditioner::BlockILU0)
            for (int p = graph.row_ptr[i]; p < diagonal_slots_[i]; ++p) {
                int j = graph.col_idx[p];
                diagonal_[j].right_solve(factors_.block(p));
                for (int update = update_offsets_[p]; update < update_offsets_[p + 1]; ++update)
                    subtract_product(factors_.block(update_targets_[update]), factors_.block(p),
                                     factors_.block(update_upper_[update]));
            }
        diagonal_[i].compute(factors_.block(diagonal_slots_[i]));
    }
}
void BlockPreconditioner::apply(const Vec &rhs, Vec &result) const {
    const auto &graph = factors_.graph;
    int width = factors_.width;
    for (int i = 0; i < graph.rows; ++i)
        work_.segment(i * width, width) = rhs.segment(order_[i] * width, width);
    for (int i = 0; i < graph.rows; ++i)
        for (int p = graph.row_ptr[i]; p < diagonal_slots_[i]; ++p)
            work_.segment(i * width, width).noalias() -=
                factors_.block(p) * work_.segment(graph.col_idx[p] * width, width);
    for (int i = graph.rows - 1; i >= 0; --i) {
        for (int p = diagonal_slots_[i] + 1; p < graph.row_ptr[i + 1]; ++p)
            work_.segment(i * width, width).noalias() -=
                factors_.block(p) * work_.segment(graph.col_idx[p] * width, width);
        diagonal_[i].solve(work_.segment(i * width, width), work_.segment(i * width, width));
    }
    for (int i = 0; i < graph.rows; ++i)
        result.segment(order_[i] * width, width) = work_.segment(i * width, width);
}
} // namespace euler

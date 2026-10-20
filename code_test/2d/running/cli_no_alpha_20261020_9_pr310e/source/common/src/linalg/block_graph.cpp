#include "euler/linalg/block_graph.hpp"
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <algorithm>
#include <numeric>
namespace euler {
BlockGraph::BlockGraph(const std::vector<std::vector<int>> &adjacency, int nx, int ny)
    : rows(static_cast<int>(adjacency.size())), cells_x(nx ? nx : rows), cells_y(ny),
      row_ptr(rows + 1, 0) {
    require(cells_x * cells_y == rows, "Block grid shape mismatch");
    for (int r = 0; r < rows; ++r) {
        auto row = adjacency[r];
        row.push_back(r);
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
        col_idx.insert(col_idx.end(), row.begin(), row.end());
        row_ptr[r + 1] = static_cast<int>(col_idx.size());
    }
}
BlockGraph BlockGraph::cartesian(int nx, int ny) {
    std::vector<std::vector<int>> adjacency(nx * ny);
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x) {
            int cell = x + nx * y;
            adjacency[cell] = {cell, (x + nx - 1) % nx + nx * y, (x + 1) % nx + nx * y,
                               x + nx * ((y + ny - 1) % ny), x + nx * ((y + 1) % ny)};
        }
    return BlockGraph(adjacency, nx, ny);
}
BlockGraph BlockGraph::periodic(int cells) {
    std::vector<std::vector<int>> adjacency(cells);
    for (int c = 0; c < cells; ++c)
        adjacency[c] = {(c + cells - 1) % cells, c, (c + 1) % cells};
    return BlockGraph(adjacency);
}
int BlockGraph::find(int row, int column) const {
    require(row >= 0 && row < rows && column >= 0 && column < rows, "Invalid block coordinates");
    auto begin = col_idx.begin() + row_ptr[row], end = col_idx.begin() + row_ptr[row + 1];
    auto found = std::lower_bound(begin, end, column);
    return found != end && *found == column ? static_cast<int>(found - col_idx.begin()) : -1;
}
BlockMatrix::BlockMatrix(const BlockGraph &g, int w, int nx_nodes, int ny_nodes)
    : graph(g), width(w), nodes_x(nx_nodes ? nx_nodes : w / NVAR), nodes_y(ny_nodes),
      values_(Vec::Zero(static_cast<Eigen::Index>(g.col_idx.size()) * w * w)) {
    require(w > 0, "Block width must be positive");
    require(nodes_x * nodes_y * NVAR == w, "Block tensor shape mismatch");
}
void BlockMatrix::zero() { values_.setZero(); }
void BlockMatrix::apply(const Vec &x, Vec &y) const {
    y.setZero();
    for (int row = 0; row < graph.rows; ++row)
        for (int p = graph.row_ptr[row]; p < graph.row_ptr[row + 1]; ++p)
            y.segment(row * width, width).noalias() +=
                block(p) * x.segment(graph.col_idx[p] * width, width);
}
std::vector<int> block_ordering(const BlockGraph &graph, bool colamd) {
    std::vector<int> order(graph.rows);
    std::iota(order.begin(), order.end(), 0);
    if (!colamd)
        return order;
    Eigen::SparseMatrix<double> connectivity(graph.rows, graph.rows);
    std::vector<Eigen::Triplet<double>> entries;
    for (int r = 0; r < graph.rows; ++r)
        for (int p = graph.row_ptr[r]; p < graph.row_ptr[r + 1]; ++p)
            entries.emplace_back(r, graph.col_idx[p], 1.0);
    connectivity.setFromTriplets(entries.begin(), entries.end());
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> permutation;
    Eigen::COLAMDOrdering<int>{}(connectivity, permutation);
    for (int i = 0; i < graph.rows; ++i)
        order[i] = permutation.indices()[i];
    return order;
}
} // namespace euler

#pragma once
#include "euler/types.hpp"
#include <vector>
namespace euler {
// CSR indices refer to cells/faces, never scalar DOFs. Diagonal is always present.
class BlockGraph {
  public:
    int rows, cells_x, cells_y;
    std::vector<int> row_ptr, col_idx;
    // Sparse outer tensor axes: [Nx_row,Ny_row,Nx_col,Ny_col].
    std::array<Eigen::Index, 4> shape() const { return {cells_x, cells_y, cells_x, cells_y}; }
    static BlockGraph periodic(int cells);
    static BlockGraph cartesian(int nx, int ny);
    explicit BlockGraph(const std::vector<std::vector<int>> &adjacency, int nx = 0, int ny = 1);
    int find(int row, int column) const;
};
class BlockMatrix {
  public:
    // Each column-major block is the 1D view of
    // [kx_row+1,ky_row+1,NVAR,kx_col+1,ky_col+1,NVAR].
    const BlockGraph graph;
    const int width, nodes_x, nodes_y;
    BlockMatrix(const BlockGraph &graph, int width, int nx_nodes = 0, int ny_nodes = 1);
    auto block(int slot) {
        return Eigen::Map<Mat>(values_.data() + static_cast<Eigen::Index>(slot) * width * width,
                               width, width);
    }
    auto block(int slot) const {
        return Eigen::Map<const Mat>(
            values_.data() + static_cast<Eigen::Index>(slot) * width * width, width, width);
    }
    auto at(int row, int col) {
        int slot = graph.find(row, col);
        require(slot >= 0, "Requested structural zero block");
        return block(slot);
    }
    auto at(int row, int col) const {
        int slot = graph.find(row, col);
        require(slot >= 0, "Requested structural zero block");
        return block(slot);
    }
    // Sparse outer coordinates select an allocated block; missing blocks are
    // structural zeros and never materialized. Matrix maps below share memory
    // with this six-axis tensor view; Eigen's first axis is always contiguous.
    auto tensor_block(int row_x, int row_y, int col_x, int col_y) {
        int index = graph.find(row_x + graph.cells_x * row_y, col_x + graph.cells_x * col_y);
        require(index >= 0 && width % NVAR == 0, "Block absent or invalid Euler block width");
        return Eigen::TensorMap<Eigen::Tensor<double, 6, Eigen::ColMajor>>(
            block(index).data(), nodes_x, nodes_y, NVAR, nodes_x, nodes_y, NVAR);
    }
    auto tensor_block(int row_x, int row_y, int col_x, int col_y) const {
        int index = graph.find(row_x + graph.cells_x * row_y, col_x + graph.cells_x * col_y);
        require(index >= 0 && width % NVAR == 0, "Block absent or invalid Euler block width");
        return Eigen::TensorMap<Eigen::Tensor<const double, 6, Eigen::ColMajor>>(
            block(index).data(), nodes_x, nodes_y, NVAR, nodes_x, nodes_y, NVAR);
    }
    void zero();
    void apply(const Vec &x, Vec &y) const;

  private:
    // One Eigen-aligned numeric pool, [block_entry,nonzero_block]. A block is
    // [row_node,row_var,col_node,col_var] in column-major order.
    Vec values_;
};
// new_to_old implements a symmetric block permutation P A P^T. COLAMD operates
// on the block connectivity graph only, during construction.
std::vector<int> block_ordering(const BlockGraph &graph, bool colamd);
} // namespace euler

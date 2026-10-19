#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <array>
#include <cmath>
#include <stdexcept>

namespace euler {
inline constexpr int NVAR = 3;
enum Component { Density, Momentum, Energy };
inline constexpr int NFACES = 2;
using Vec = Eigen::VectorXd;
using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using BlockRef = Eigen::Ref<Mat>;
using ConstBlockRef = Eigen::Ref<const Mat>;
using State = Eigen::Matrix<double, NVAR, 1, Eigen::ColMajor>;
using FluxJacobian = Eigen::Matrix<double, NVAR, NVAR, Eigen::ColMajor>;
inline void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
// Column-major tensor [kx+1,ky+1,NVAR,Nx,Ny], first axis contiguous.
// In 1D: [nodes,1,NVAR,cells,1]. Each cell is a zero-copy Eigen
// column-major [nodes,NVAR] map; node varies fastest, then variable, then cell.
struct Layout {
    int cells, nodes;
    std::array<Eigen::Index, 5> shape() const { return {nodes, 1, NVAR, cells, 1}; }
    auto tensor(Vec &u) const {
        return Eigen::TensorMap<Eigen::Tensor<double, 5, Eigen::ColMajor>>(u.data(), nodes, 1, NVAR,
                                                                           cells, 1);
    }
    auto tensor(const Vec &u) const {
        return Eigen::TensorMap<Eigen::Tensor<const double, 5, Eigen::ColMajor>>(u.data(), nodes, 1,
                                                                                 NVAR, cells, 1);
    }
    int block_size() const { return NVAR * nodes; }
    int size() const { return cells * block_size(); }
    int local(int variable, int node) const { return variable * nodes + node; }
    int index(int cell, int variable, int node) const {
        return cell * block_size() + local(variable, node);
    }
    auto cell(Vec &u, int c) const {
        return Eigen::Map<Mat>(u.data() + c * block_size(), nodes, NVAR);
    }
    auto cell(const Vec &u, int c) const {
        return Eigen::Map<const Mat>(u.data() + c * block_size(), nodes, NVAR);
    }
};
} // namespace euler

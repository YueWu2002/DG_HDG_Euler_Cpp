#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <array>
#include <cmath>
#include <stdexcept>
namespace euler {
inline constexpr int NDIM = 2;
inline constexpr int NVAR = NDIM + 2;
inline constexpr int NFACES = 2 * NDIM;
enum Component { Density, MomentumX, MomentumY, Energy };
using Vec = Eigen::VectorXd;
using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using BlockRef = Eigen::Ref<Mat>;
using ConstBlockRef = Eigen::Ref<const Mat>;
using State = Eigen::Matrix<double, NVAR, 1>;
using FluxJacobian = Eigen::Matrix<double, NVAR, NVAR>;
inline void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
struct Layout {
    int nx, ny, nodes, cells;
    Layout(int nx_, int ny_, int n) : nx(nx_), ny(ny_), nodes(n), cells(nx_ * ny_) {}
    int plane() const { return nodes * nodes; }
    int block_size() const { return NVAR * plane(); }
    int size() const { return cells * block_size(); }
    int local(int variable, int node) const { return variable * plane() + node; }
    int index(int cell, int variable, int node) const {
        return cell * block_size() + local(variable, node);
    }
    std::array<Eigen::Index, 5> shape() const { return {nodes, nodes, NVAR, nx, ny}; }
    auto tensor(Vec &u) const {
        return Eigen::TensorMap<Eigen::Tensor<double, 5, Eigen::ColMajor>>(
            u.data(), nodes, nodes, NVAR, nx, ny);
    }
    auto tensor(const Vec &u) const {
        return Eigen::TensorMap<Eigen::Tensor<const double, 5, Eigen::ColMajor>>(
            u.data(), nodes, nodes, NVAR, nx, ny);
    }
    auto cell(Vec &u, int c) const {
        return Eigen::Map<Mat>(u.data() + c * block_size(), plane(), NVAR);
    }
    auto cell(const Vec &u, int c) const {
        return Eigen::Map<const Mat>(u.data() + c * block_size(), plane(), NVAR);
    }
    auto sheet(Vec &u, int c, int variable) const {
        return Eigen::Map<Mat>(u.data() + index(c, variable, 0), nodes, nodes);
    }
    auto sheet(const Vec &u, int c, int variable) const {
        return Eigen::Map<const Mat>(u.data() + index(c, variable, 0), nodes, nodes);
    }
};
} // namespace euler

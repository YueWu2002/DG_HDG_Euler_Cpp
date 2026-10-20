#pragma once
#include "euler/types.hpp"
#include <limits>
namespace euler {
// Continuous LU pool and one factorization workspace. No duplicate BCSR
// diagonals, and no per-block dense matrices or triangular-solve allocations.
class BlockJacobi {
  public:
    BlockJacobi(int blocks,int width)
        : blocks_(blocks),width_(width),factors_(static_cast<Eigen::Index>(blocks)*width*width),
          permutations_(blocks*width),block_(width,width),lu_(width),work_(width) {}
    template<class Provider> void refresh(Provider &&provider) {
        for(int i=0;i<blocks_;++i) {
            provider(i,i,block_);
            lu_.compute(block_);
            require(lu_.matrixLU().allFinite() &&
                (lu_.matrixLU().diagonal().array().abs()>std::numeric_limits<double>::min()).all(),
                "Singular block Jacobi factor");
            Eigen::Map<Mat> factor(factors_.data()+static_cast<Eigen::Index>(i)*width_*width_,width_,width_);
            factor=lu_.matrixLU();
            permutations_.segment(i*width_,width_)=lu_.permutationP().indices();
        }
    }
    void apply(const Vec &rhs,Vec &result) const {
        for(int i=0;i<blocks_;++i) {
            const Eigen::Map<const Mat> factor(factors_.data()+static_cast<Eigen::Index>(i)*width_*width_,width_,width_);
            for(int j=0;j<width_;++j) work_[permutations_[i*width_+j]]=rhs[i*width_+j];
            if(width_==NVAR) {
                Eigen::Map<const FluxJacobian> small_factor(factor.data());
                Eigen::Map<State> small_rhs(work_.data());
                small_factor.template triangularView<Eigen::UnitLower>().solveInPlace(small_rhs);
                small_factor.template triangularView<Eigen::Upper>().solveInPlace(small_rhs);
            } else {
                factor.template triangularView<Eigen::UnitLower>().solveInPlace(work_);
                factor.template triangularView<Eigen::Upper>().solveInPlace(work_);
            }
            result.segment(i*width_,width_)=work_;
        }
    }
  private:
    int blocks_,width_;
    Vec factors_;
    Eigen::VectorXi permutations_;
    Mat block_;
    Eigen::PartialPivLU<Mat> lu_;
    mutable Vec work_;
};
}

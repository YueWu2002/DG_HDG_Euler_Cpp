#include "euler/linalg/gmres.hpp"
#include <algorithm>
#include <limits>
namespace euler {
GMRES::GMRES(int n, int restart, const Vec &inverse_mass)
    : restart_(std::min(n, restart)), basis_(n, restart_ + 1),
      hessenberg_(restart_ + 1, restart_), cosine_(restart_), sine_(restart_),
      rotated_rhs_(restart_ + 1), coefficients_(restart_), work_(n), residual_(n), column_(n),
      weights_(inverse_mass), weighted_(n), projection_(restart_) {
    require(n > 0 && restart > 0 && weights_.size() == n && weights_.allFinite() &&
            (weights_.array() > 0).all(), "Invalid GMRES metric or dimensions");
}
double GMRES::dot(const Vec &a, const Vec &b) const {
    return (a.array() * weights_.array() * b.array()).sum();
}
double GMRES::norm(const Vec &a) const { return std::sqrt(std::max(0.0, dot(a, a))); }
LinearResult GMRES::solve(Action A, Action M, const Vec &rhs, Vec &x, int limit, double tolerance) {
    x.setZero();
    residual_ = rhs;
    double rhs_norm = norm(rhs);
    if (rhs_norm == 0)
        return {0, 0};
    int iterations = 0;
    while (iterations < limit) {
        double beta = norm(residual_);
        if (beta <= tolerance * rhs_norm)
            return {iterations, beta / rhs_norm};
        basis_.col(0) = residual_ / beta;
        hessenberg_.setZero();
        rotated_rhs_.setZero();
        rotated_rhs_[0] = beta;
        int count = 0;
        for (int j = 0; j < restart_ && iterations < limit; ++j) {
            column_ = basis_.col(j);
            M(column_, work_);
            A(work_, column_);
            // CGS2: two weighted projections and two fused basis actions.
            // The fixed right preconditioner is applied to the combined update.
            for (int pass = 0; pass < 2; ++pass) {
                weighted_.array() = weights_.array() * column_.array();
                projection_.head(j + 1).noalias() = basis_.leftCols(j + 1).transpose() * weighted_;
                hessenberg_.col(j).head(j + 1) += projection_.head(j + 1);
                work_.noalias() = basis_.leftCols(j + 1) * projection_.head(j + 1);
                column_ -= work_;
            }
            double next = norm(column_);
            hessenberg_(j + 1, j) = next;
            if (next > 0)
                basis_.col(j + 1) = column_ / next;
            for (int i = 0; i < j; ++i) {
                double a = hessenberg_(i, j), b = hessenberg_(i + 1, j);
                hessenberg_(i, j) = cosine_[i] * a + sine_[i] * b;
                hessenberg_(i + 1, j) = -sine_[i] * a + cosine_[i] * b;
            }
            double d = std::hypot(hessenberg_(j, j), hessenberg_(j + 1, j));
            require(std::isfinite(d) && d > 0, "GMRES Arnoldi breakdown");
            cosine_[j] = hessenberg_(j, j) / d;
            sine_[j] = hessenberg_(j + 1, j) / d;
            hessenberg_(j, j) = d;
            hessenberg_(j + 1, j) = 0;
            rotated_rhs_[j + 1] = -sine_[j] * rotated_rhs_[j];
            rotated_rhs_[j] *= cosine_[j];
            ++iterations;
            count = j + 1;
            if (std::abs(rotated_rhs_[j + 1]) <= tolerance * rhs_norm ||
                next <= std::numeric_limits<double>::min())
                break;
        }
        for (int i = count - 1; i >= 0; --i) {
            double value = rotated_rhs_[i];
            for (int j = i + 1; j < count; ++j)
                value -= hessenberg_(i, j) * coefficients_[j];
            coefficients_[i] = value / hessenberg_(i, i);
        }
        work_.noalias() = basis_.leftCols(count) * coefficients_.head(count);
        M(work_, column_);
        x += column_;
        A(x, work_);
        residual_ = rhs - work_;
        double relative = norm(residual_) / rhs_norm;
        require(std::isfinite(relative), "Nonfinite GMRES residual");
        if (relative <= tolerance)
            return {iterations, relative};
    }
    throw std::runtime_error("GMRES iteration limit reached");
}
} // namespace euler

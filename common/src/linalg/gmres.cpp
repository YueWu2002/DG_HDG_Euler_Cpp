#include "euler/linalg/gmres.hpp"
#include <algorithm>
#include <limits>
namespace euler {
GMRES::GMRES(int n, int restart, const Vec &inverse_mass)
    : restart_(std::min(n, restart)), basis_(n, restart_ + 1), preconditioned_(n, restart_),
      hessenberg_(restart_ + 1, restart_), cosine_(restart_), sine_(restart_),
      rotated_rhs_(restart_ + 1), coefficients_(restart_), work_(n), residual_(n), column_(n),
      weights_(inverse_mass) {}
double GMRES::dot(const Vec &a, const Vec &b) const {
    return (a.array() * weights_.array() * b.array()).sum();
}
double GMRES::norm(const Vec &a) const { return std::sqrt(std::max(0.0, dot(a, a))); }
LinearResult GMRES::solve(Action A, Action M, const Vec &rhs, Vec &x, int limit, double tolerance) {
    x.setZero();
    double rhs_norm = norm(rhs);
    if (rhs_norm == 0)
        return {0, 0};
    residual_ = rhs;
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
            preconditioned_.col(j) = work_;
            A(work_, column_);
            for (int pass = 0; pass < 2; ++pass)
                for (int i = 0; i <= j; ++i) {
                    work_ = basis_.col(i);
                    double h = dot(work_, column_);
                    hessenberg_(i, j) += h;
                    column_.noalias() -= h * work_;
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
        for (int i = 0; i < count; ++i)
            x.noalias() += coefficients_[i] * preconditioned_.col(i);
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

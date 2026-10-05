#include "residual.h"
#include <cmath>
#include <vector>

static double frobenius_norm(const double* M, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n * n; ++i) s += M[i] * M[i];
    return std::sqrt(s);
}

double inverse_residual(const double* A, const double* X, std::size_t n) {
    std::vector<double> R(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                s += A[i * n + k] * X[k * n + j];
            }
            R[i * n + j] = s - ((i == j) ? 1.0 : 0.0);
        }
    }
    return frobenius_norm(R.data(), n);
}
#include "jordan.h"
#include <cmath>
#include <limits>
#include <vector>
#include <utility>

bool jordan_inverse_col_pivot(double* A, std::size_t n) {
    std::vector<double> E(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) E[i * n + i] = 1.0;

    const double tiny = std::numeric_limits<double>::min();

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot_row = k;
        double max_abs = std::fabs(A[k * n + k]);
        for (std::size_t i = k + 1; i < n; ++i) {
            double v = std::fabs(A[i * n + k]);
            if (v > max_abs) { max_abs = v; pivot_row = i; }
        }
        if (max_abs < tiny) return false;

        if (pivot_row != k) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(A[k * n + j], A[pivot_row * n + j]);
                std::swap(E[k * n + j], E[pivot_row * n + j]);
            }
        }

        double inv_diag = 1.0 / A[k * n + k];
        for (std::size_t j = 0; j < n; ++j) {
            A[k * n + j] *= inv_diag;
            E[k * n + j] *= inv_diag;
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            double factor = A[i * n + k];
            if (factor == 0.0) continue;
            for (std::size_t j = 0; j < n; ++j) {
                A[i * n + j] -= factor * A[k * n + j];
                E[i * n + j] -= factor * E[k * n + j];
            }
        }
    }

    for (std::size_t i = 0; i < n * n; ++i) A[i] = E[i];
    return true;
}
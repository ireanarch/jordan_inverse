#include "matrix_init.h"
#include <cmath>

bool init_matrix_formula(double* A, std::size_t n, int k) {
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double val = 0.0;
            switch (k) {
            case 1: {
                std::size_t mx = (i > j) ? i : j;
                val = static_cast<double>(n - mx);
                break;
            }
            case 2: {
                std::size_t mx = (i > j) ? i : j;
                val = static_cast<double>(mx + 1);
                break;
            }
            case 3:
                val = std::fabs(static_cast<double>(i) - static_cast<double>(j));
                break;
            case 4:
                val = 1.0 / static_cast<double>(i + j + 1);
                break;
            default:
                return false;
            }
            A[i * n + j] = val;
        }
    }
    return true;
}
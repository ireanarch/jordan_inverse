#ifndef RESIDUAL_H
#define RESIDUAL_H

#include <cstddef>

// Последовательная версия
double inverse_residual(const double* A, const double* X, std::size_t n);

// Параллельная версия
double inverse_residual_pthreads(const double* A, const double* X, std::size_t n, int num_threads);

#endif // RESIDUAL_H
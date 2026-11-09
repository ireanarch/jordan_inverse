#ifndef JORDAN_PTHREADS_H
#define JORDAN_PTHREADS_H

#include <cstddef>

// Возвращает true при успехе, false если матрица вырождена
bool jordan_inverse_pthreads(double* A, double* A_inv, std::size_t n, int num_threads);

#endif // JORDAN_PTHREADS_H

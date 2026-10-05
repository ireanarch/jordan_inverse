#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include <cstddef>
#include <string>

bool read_matrix_file(const std::string& filename, double* A, std::size_t n);

void print_matrix(const double* M, std::size_t l, std::size_t n, std::size_t m);

#endif // MATRIX_IO_H
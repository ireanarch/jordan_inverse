#include "matrix_io.h"
#include <cstdio>
#include <fstream>

bool read_matrix_file(const std::string& filename, double* A, std::size_t n) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::fprintf(stderr, "Error: cannot open file '%s'\n", filename.c_str());
        return false;
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (!(fin >> A[i * n + j])) {
                std::fprintf(stderr,
                    "Error: not enough data or invalid format in file '%s'\n",
                    filename.c_str());
                return false;
            }
        }
    }
    return true;
}

void print_matrix(const double* M, std::size_t l, std::size_t n, std::size_t m) {
    std::size_t rows = (l < m) ? l : m;
    std::size_t cols = (n < m) ? n : m;
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            std::printf(" %10.3e", M[i * n + j]);
        }
        std::printf("\n");
    }
}
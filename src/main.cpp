#include "matrix_init.h"
#include "matrix_io.h"
#include "jordan.h"
#include "residual.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::fprintf(stderr,
            "Usage: %s n m k [filename]\n"
            "  n        - matrix dimension\n"
            "  m        - max rows/cols to print\n"
            "  k        - formula number (1..4) or 0 (read from file)\n"
            "  filename - required only if k == 0\n",
            argv[0]);
        return 1;
    }

    std::size_t n = static_cast<std::size_t>(std::strtoull(argv[1], nullptr, 10));
    std::size_t m = static_cast<std::size_t>(std::strtoull(argv[2], nullptr, 10));
    int k = std::atoi(argv[3]);

    if (n == 0 || m == 0) {
        std::fprintf(stderr, "Error: n and m must be positive\n");
        return 1;
    }
    if (k < 0 || k > 4) {
        std::fprintf(stderr, "Error: k must be 0, 1, 2, 3 or 4\n");
        return 1;
    }

    std::string filename;
    if (k == 0) {
        if (argc < 5) {
            std::fprintf(stderr, "Error: filename required when k == 0\n");
            return 1;
        }
        filename = argv[4];
    }

    std::vector<double> A(n * n);

    if (k == 0) {
        if (!read_matrix_file(filename, A.data(), n)) return 1;
    }
    else {
        if (!init_matrix_formula(A.data(), n, k)) {
            std::fprintf(stderr, "Error: unknown formula number %d\n", k);
            return 1;
        }
    }

    std::printf("\nInitial matrix A:\n");
    print_matrix(A.data(), n, n, m);

    std::clock_t t0 = std::clock();
    bool ok = jordan_inverse_col_pivot(A.data(), n);
    std::clock_t t1 = std::clock();

    if (!ok) {
        std::fprintf(stderr, "Error: matrix is singular\n");
        return 1;
    }

    double elapsed = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
    std::printf("\nInverse matrix X:\n");
    print_matrix(A.data(), n, n, m);

    std::vector<double> A2(n * n);
    if (k == 0) {
        if (!read_matrix_file(filename, A2.data(), n)) return 1;
    }
    else {
        if (!init_matrix_formula(A2.data(), n, k)) return 1;
    }

    double res = inverse_residual(A2.data(), A.data(), n);
    std::printf("\nResidual norm ||A*A^{-1} - E|| = %.3e\n", res);
    std::printf("Elapsed time: %.2f s\n", elapsed);

    return 0;
}
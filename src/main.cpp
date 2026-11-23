#include "matrix_init.h"
#include "matrix_io.h"
#include "jordan.h"
#include "jordan_pthreads.h"
#include "residual.h"

#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

static double elapsed_seconds(
    const std::chrono::steady_clock::time_point& t0,
    const std::chrono::steady_clock::time_point& t1)
{
    return std::chrono::duration<double>(t1 - t0).count();
}

int main(int argc, char* argv[]) {
    auto total_start = std::chrono::steady_clock::now();

    if (argc < 4) {
        std::fprintf(stderr,
            "Usage: %s n m k [filename] [num_threads]\n"
            "  n           - matrix dimension\n"
            "  m           - max rows/cols to print\n"
            "  k           - formula number (1..4) or 0 (read from file)\n"
            "  filename    - required only if k == 0\n"
            "  num_threads - number of threads (default: 1)\n",
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

    // ==========================================
    // ЧТЕНИЕ ЧИСЛА ПОТОКОВ
    // ==========================================
    int num_threads = 1; 
    if (k == 0) {
        // Если файл, то потоки идут 6-м аргументом (argc >= 6)
        if (argc >= 6) num_threads = std::atoi(argv[5]);
    } else {
        // Если формула, то потоки идут 5-м аргументом (argc >= 5)
        if (argc >= 5) num_threads = std::atoi(argv[4]);
    }
    
    if (num_threads < 1) num_threads = 1;
    
    std::printf("Running with %d thread(s)...\n", num_threads);

    // ==========================================
    // ВЫДЕЛЕНИЕ ПАМЯТИ И ИНИЦИАЛИЗАЦИЯ
    // ==========================================
    auto t_alloc_start = std::chrono::steady_clock::now();
    std::vector<double> A(n * n);
    std::vector<double> A_inv(n * n); // Матрица для результата
    std::vector<double> A2(n * n);    // Для повторной инициализации (невязка)
    auto t_alloc_end = std::chrono::steady_clock::now();
    double time_alloc = elapsed_seconds(t_alloc_start, t_alloc_end);

    auto t_init_start = std::chrono::steady_clock::now();
    if (k == 0) {
        if (!read_matrix_file(filename, A.data(), n)) return 1;
    } else {
        if (!init_matrix_formula(A.data(), n, k)) {
            std::fprintf(stderr, "Error: unknown formula number %d\n", k);
            return 1;
        }
    }
    auto t_init_end = std::chrono::steady_clock::now();
    double time_init = elapsed_seconds(t_init_start, t_init_end);

    auto t_print1_start = std::chrono::steady_clock::now();
    std::printf("\nInitial matrix A:\n");
    print_matrix(A.data(), n, n, m);
    auto t_print1_end = std::chrono::steady_clock::now();
    double time_print1 = elapsed_seconds(t_print1_start, t_print1_end);

    // ==========================================
    // НАХОЖДЕНИЕ ОБРАТНОЙ МАТРИЦЫ
    // ==========================================
    auto t_jordan_start = std::chrono::steady_clock::now();
    bool ok = false;
    
    if (num_threads > 1) {
        // Вызов параллельной версии
        ok = jordan_inverse_pthreads(A.data(), A_inv.data(), n, num_threads);
    } else {
        // Вызов последовательной версии
        ok = jordan_inverse_col_pivot(A.data(), n);
        // Так как старая функция записывает результат в A, копируем его в A_inv
        std::copy(A.begin(), A.end(), A_inv.begin());
    }
    
    auto t_jordan_end = std::chrono::steady_clock::now();
    double time_jordan = elapsed_seconds(t_jordan_start, t_jordan_end);

    if (!ok) {
        std::fprintf(stderr, "Error: matrix is singular\n");
        return 1;
    }

    auto t_print2_start = std::chrono::steady_clock::now();
    std::printf("\nInverse matrix X:\n");
    print_matrix(A_inv.data(), n, n, m);
    auto t_print2_end = std::chrono::steady_clock::now();
    double time_print2 = elapsed_seconds(t_print2_start, t_print2_end);

    // ==========================================
    // ПОВТОРНАЯ ИНИЦИАЛИЗАЦИЯ И НЕВЯЗКА
    // ==========================================
    auto t_reinit_start = std::chrono::steady_clock::now();
    if (k == 0) {
        if (!read_matrix_file(filename, A2.data(), n)) return 1;
    } else {
        if (!init_matrix_formula(A2.data(), n, k)) return 1;
    }
    auto t_reinit_end = std::chrono::steady_clock::now();
    double time_reinit = elapsed_seconds(t_reinit_start, t_reinit_end);

    auto t_res_start = std::chrono::steady_clock::now();
    double res;
    if (num_threads > 1) {
        res = inverse_residual_pthreads(A2.data(), A_inv.data(), n, num_threads);
    }
    else {
        res = inverse_residual(A2.data(), A_inv.data(), n);
    }
    auto t_res_end = std::chrono::steady_clock::now();
    double time_res = elapsed_seconds(t_res_start, t_res_end);

    auto total_end = std::chrono::steady_clock::now();
    double time_total = elapsed_seconds(total_start, total_end);

    // ==========================================
    // ВЫВОД ТАЙМИНГОВ
    // ==========================================
    std::printf("\nResidual norm ||A*A^{-1} - E|| = %.3e\n\n", res);

    std::printf("========== Timings ==========\n");
    std::printf("Allocation             : %.3f s\n", time_alloc);
    std::printf("Init matrix A          : %.3f s\n", time_init);
    std::printf("Print initial matrix   : %.3f s\n", time_print1);
    std::printf("Jordan inverse         : %.3f s\n", time_jordan);
    std::printf("Print inverse matrix   : %.3f s\n", time_print2);
    std::printf("Re-init matrix A       : %.3f s\n", time_reinit);
    std::printf("Residual norm          : %.3f s\n", time_res);
    std::printf("--------------------------------\n");
    std::printf("TOTAL                  : %.3f s\n", time_total);
    std::printf("=============================\n");
    return 0;
}

#include "residual.h"
#include <cmath>
#include <vector>
#include <pthread.h>

// Вспомогательная функция для нормы Фробениуса
static double frobenius_norm(const double* M, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n * n; ++i) s += M[i] * M[i];
    return std::sqrt(s);
}

// Старая последовательная версия
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

// Структура для передачи данных в поток
struct ResThreadData {
    const double* A;
    const double* X;
    double* R;
    std::size_t n;
    std::size_t start_row;
    std::size_t end_row;
};

// Функция, которую выполняет каждый поток
void* residual_worker(void* arg) {
    ResThreadData* data = (ResThreadData*)arg;
    std::size_t n = data->n;

    for (std::size_t i = data->start_row; i < data->end_row; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                s += data->A[i * n + k] * data->X[k * n + j];
            }
            data->R[i * n + j] = s - ((i == j) ? 1.0 : 0.0);
        }
    }
    return nullptr;
}

// Новая параллельная версия
double inverse_residual_pthreads(const double* A, const double* X, std::size_t n, int num_threads) {
    if (num_threads <= 0) num_threads = 1;
    if (num_threads > (int)n) num_threads = (int)n; // Больше потоков, чем строк, не нужно

    std::vector<double> R(n * n, 0.0);
    std::vector<pthread_t> threads(num_threads);
    std::vector<ResThreadData> thread_data(num_threads);

    std::size_t rows_per_thread = n / num_threads;

    for (int t = 0; t < num_threads; ++t) {
        std::size_t start_row = t * rows_per_thread;
        std::size_t end_row = (t == num_threads - 1) ? n : start_row + rows_per_thread;

        thread_data[t] = { A, X, R.data(), n, start_row, end_row };

        int rc = pthread_create(&threads[t], nullptr, residual_worker, &thread_data[t]);
        if (rc != 0) {
            // Если поток не создался — считаем последовательно (страховка)
            return inverse_residual(A, X, n);
        }
    }

    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], nullptr);
    }

    return frobenius_norm(R.data(), n);
}
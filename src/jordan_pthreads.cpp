#include "jordan_pthreads.h"
#include <cmath>
#include <limits>
#include <vector>
#include <pthread.h>
#include <iostream>

struct ThreadData {
    double* A;
    std::size_t n;
    int thread_id;
    int num_threads;
    int k;
    pthread_barrier_t* barrier;
};

void* jordan_worker(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    std::size_t n = data->n;
    double* aug = data->A; 
    int p = data->num_threads;
    int tid = data->thread_id;
    const double tiny = std::numeric_limits<double>::min();
    
    for (std::size_t k = 0; k < n; ++k) {
        data->k = k;
        
        // ЭТАП 1: Подготовка (только 0-й поток)
        if (tid == 0) {
            std::size_t pivot_row = k;
            double max_abs = std::fabs(aug[k * 2 * n + k]);
            for (std::size_t i = k + 1; i < n; ++i) {
                double v = std::fabs(aug[i * 2 * n + k]);
                if (v > max_abs) { max_abs = v; pivot_row = i; }
            }
            
            if (max_abs < tiny) {
                aug[0] = std::numeric_limits<double>::quiet_NaN(); // Флаг вырожденности
            } else {
                if (pivot_row != k) {
                    for (std::size_t j = 0; j < 2 * n; ++j) {
                        std::swap(aug[k * 2 * n + j], aug[pivot_row * 2 * n + j]);
                    }
                }
                double inv_diag = 1.0 / aug[k * 2 * n + k];
                for (std::size_t j = 0; j < 2 * n; ++j) {
                    aug[k * 2 * n + j] *= inv_diag;
                }
            }
        }
        
        pthread_barrier_wait(data->barrier); // Ждем подготовку
        
        if (std::isnan(aug[0])) break;
        
        // ЭТАП 2: Исключение (все потоки)
        for (std::size_t i = tid; i < n; i += p) {
            if (i == k) continue;
            double factor = aug[i * 2 * n + k];
            if (std::fabs(factor) < tiny) continue;
            for (std::size_t j = 0; j < 2 * n; ++j) {
                aug[i * 2 * n + j] -= factor * aug[k * 2 * n + j];
            }
        }
        
        pthread_barrier_wait(data->barrier); // Ждем окончания исключения
    }
    return nullptr;
}

bool jordan_inverse_pthreads(double* A, double* A_inv, std::size_t n, int num_threads) {
    if (num_threads <= 0) num_threads = 1;
    
    std::vector<double> aug(n * 2 * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            aug[i * 2 * n + j] = A[i * n + j];
            aug[i * 2 * n + (n + j)] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, nullptr, num_threads);
    
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> thread_data(num_threads);
    
    for (int t = 0; t < num_threads; ++t) {
        thread_data[t] = { aug.data(), n, t, num_threads, 0, &barrier };
        pthread_create(&threads[t], nullptr, jordan_worker, &thread_data[t]);
    }
    
    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], nullptr);
    }
    
    pthread_barrier_destroy(&barrier);
    
    if (std::isnan(aug[0])) return false;
    
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            A_inv[i * n + j] = aug[i * 2 * n + (n + j)];
        }
    }
    return true;
}

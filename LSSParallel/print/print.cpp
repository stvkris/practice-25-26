
#include <cstdlib>
#include <pthread.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <limits>

struct ThreadData {
    double* A;
    double* x;
    double* b;
    int n;
    int start_row;
    int end_row;
    double partial_norm_top;
    double partial_norm_bot;
};

void print_m(double* A, int n, int m)
{
    int s = std::min(n, m);
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            std::cout << A[i * n + j] << " ";
        }
        std::cout << "\n";
    }
}

void print_v(double* b, int n, int m)
{
    int s = std::min(n, m);
    for (int i = 0; i < s; i++) {
        std::cout << b[i] << "\n";
    }
}

void* compute_partial_norm(void* arg) {
    ThreadData* data = static_cast<ThreadData*>(arg);
    
    data->partial_norm_top = 0.0;
    data->partial_norm_bot = 0.0;
    
    for (int i = data->start_row; i < data->end_row; i++) {
        double sum = 0.0;
        for (int j = 0; j < data->n; j++) {
            sum += data->A[i * data->n + j] * data->x[j];
        }
        double diff = sum - data->b[i];
        data->partial_norm_top += diff * diff;
        data->partial_norm_bot += data->b[i] * data->b[i];
    }
    
    return nullptr;
}

double norm(double* A, double* x, double* b, int n, int num_threads)
{
    double epsilon_d = std::numeric_limits<double>::epsilon();
    
    if (n < num_threads) {
        num_threads = n;
    }
    
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> thread_data(num_threads);
        
    int rows_per_thread = n / num_threads;
    int remaining_rows = n % num_threads;
    
    int current_row = 0;
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].A = A;
        thread_data[i].x = x;
        thread_data[i].b = b;
        thread_data[i].n = n;
        thread_data[i].start_row = current_row;
        
        int rows_for_this_thread = rows_per_thread;
        if (i < remaining_rows) {
            rows_for_this_thread++;
        }
        
        thread_data[i].end_row = current_row + rows_for_this_thread;
        current_row = thread_data[i].end_row;
        
        pthread_create(&threads[i], nullptr, compute_partial_norm, &thread_data[i]);
    }
    
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], nullptr);
    }

    double norm_top = 0.0;
    double norm_bot = 0.0;
    
    for (int i = 0; i < num_threads; i++) {
        norm_top += thread_data[i].partial_norm_top;
        norm_bot += thread_data[i].partial_norm_bot;
    }
    
    if (norm_bot <= epsilon_d * n) { 
        return -1.0;
    }
    return sqrt(norm_top / norm_bot);
}


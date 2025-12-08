#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <pthread.h>
#include <iostream>
#include <queue>
#include <functional>
#include <stdexcept>

struct ThreadData {
    double* A;
    double* b;
    int n;
    int start_row;
    int end_row;
    int current_pivot;
};

struct PivotThreadData {
    double* A;
    int n;
    int start_row;
    int end_row;
    int current_pivot;
    int local_pivot_row;
    double local_max_val;
};

struct BackSubThreadData {
    double* A;
    double* x;
    int n;
    int i;
    int start_j;
    int end_j;
    double* partial_sum;
};

class ThreadPool {
private:
    std::vector<pthread_t> threads;
    std::queue<std::function<void()>> tasks;
    pthread_mutex_t queue_mutex;
    pthread_cond_t condition;
    bool stop;
    int task_counter;
    pthread_mutex_t counter_mutex;
    pthread_cond_t all_tasks_done;
    
    static void* worker(void* arg) {
        ThreadPool* pool = static_cast<ThreadPool*>(arg);
        pool->worker_thread();
        return nullptr;
    }
    
    void worker_thread() {
        while (true) {
            std::function<void()> task;
            
            pthread_mutex_lock(&queue_mutex);
            while (tasks.empty() && !stop) {
                pthread_cond_wait(&condition, &queue_mutex);
            }
            
            if (stop && tasks.empty()) {
                pthread_mutex_unlock(&queue_mutex);
                return;
            }
            
            task = std::move(tasks.front());
            tasks.pop();
            pthread_mutex_unlock(&queue_mutex);
            
            task();
            
            pthread_mutex_lock(&counter_mutex);
            task_counter--;
            if (task_counter == 0) {
                pthread_cond_signal(&all_tasks_done);
            }
            pthread_mutex_unlock(&counter_mutex);
        }
    }
    
public:
    ThreadPool(size_t num_threads) : stop(false), task_counter(0) {
        pthread_mutex_init(&queue_mutex, nullptr);
        pthread_cond_init(&condition, nullptr);
        pthread_mutex_init(&counter_mutex, nullptr);
        pthread_cond_init(&all_tasks_done, nullptr);
        
        threads.resize(num_threads);
        for (size_t i = 0; i < num_threads; ++i) {
            pthread_create(&threads[i], nullptr, &ThreadPool::worker, this);
        }
    }
    
    ~ThreadPool() {
        pthread_mutex_lock(&queue_mutex);
        stop = true;
        pthread_mutex_unlock(&queue_mutex);
        
        pthread_cond_broadcast(&condition);
        
        for (pthread_t& thread : threads) {
            pthread_join(thread, nullptr);
        }
        
        pthread_mutex_destroy(&queue_mutex);
        pthread_cond_destroy(&condition);
        pthread_mutex_destroy(&counter_mutex);
        pthread_cond_destroy(&all_tasks_done);
    }
    
    void wait_all() {
        pthread_mutex_lock(&counter_mutex);
        while (task_counter > 0) {
            pthread_cond_wait(&all_tasks_done, &counter_mutex);
        }
        pthread_mutex_unlock(&counter_mutex);
    }
    
    template<class F>
    void enqueue(F&& task) {
        pthread_mutex_lock(&counter_mutex);
        task_counter++;
        pthread_mutex_unlock(&counter_mutex);
        
        pthread_mutex_lock(&queue_mutex);
        tasks.push(std::forward<F>(task));
        pthread_mutex_unlock(&queue_mutex);
        pthread_cond_signal(&condition);
    }
    
    size_t size() const {
        return threads.size();
    }
};

static ThreadPool* global_thread_pool = nullptr;

void init_thread_pool(int num_threads) {
    if (global_thread_pool == nullptr) {
        global_thread_pool = new ThreadPool(num_threads);
    }
}

void cleanup_thread_pool() {
    if (global_thread_pool != nullptr) {
        delete global_thread_pool;
        global_thread_pool = nullptr;
    }
}

void find_pivot_row_partial_task(PivotThreadData* data) {
    data->local_max_val = fabs(data->A[data->current_pivot * data->n + data->current_pivot]);
    data->local_pivot_row = data->current_pivot;
    
    for (int k = data->start_row; k < data->end_row; k++) {
        if (k == data->current_pivot) continue;
        
        double val = fabs(data->A[k * data->n + data->current_pivot]);
        if (val > data->local_max_val) {
            data->local_max_val = val;
            data->local_pivot_row = k;
        }
    }
}

void null_rows_task(ThreadData* data) {
    for (int k = data->start_row; k < data->end_row; k++) {
        if (k == data->current_pivot) continue;
        
        double factor = data->A[k * data->n + data->current_pivot] / data->A[data->current_pivot * data->n + data->current_pivot];
        for (int j = data->current_pivot; j < data->n; j++) {
            data->A[k * data->n + j] -= data->A[data->current_pivot * data->n + j] * factor;
        }
        data->b[k] -= data->b[data->current_pivot] * factor;
    }
}

void compute_partial_sum_task(BackSubThreadData* data) {
    double sum = 0.0;
    for (int j = data->start_j; j < data->end_j; j++) {
        sum += data->A[data->i * data->n + j] * data->x[j];
    }
    *data->partial_sum = sum;
}

bool solve(double* A, double* b, int n, double* x, int num_threads)
{
    if (n <= 0 || num_threads <= 0) {
        return false;
    }
    
    if (global_thread_pool == nullptr) {
        init_thread_pool(num_threads);
    }
    
    num_threads = global_thread_pool->size();
    if (n < num_threads) {
        num_threads = n;
    }
    
    double epsilon_d = std::numeric_limits<double>::epsilon();
    
    double* A_copy = new double[n * n];
    double* b_copy = new double[n];
    
    for (int i = 0; i < n * n; i++) {
        A_copy[i] = A[i];
    }
    for (int i = 0; i < n; i++) {
        b_copy[i] = b[i];
        x[i] = 0.0;
    }
    
    bool use_multithreading = (n > 100);
    
    for (int i = 0; i < n; i++) {
        int pivotRow = i;
        double maxVal = fabs(A_copy[i * n + i]);
        
        if (n - i > 1) {
            if (use_multithreading && (n - i - 1 > 100)) {
                std::vector<PivotThreadData> thread_data(num_threads);
                
                int rows_per_thread = (n - i) / num_threads;
                int remaining_rows = (n - i) % num_threads;
                
                int current_row = i;
                for (int t = 0; t < num_threads; t++) {
                    thread_data[t].A = A_copy;
                    thread_data[t].n = n;
                    thread_data[t].current_pivot = i;
                    thread_data[t].start_row = current_row;
                    
                    int rows_for_this_thread = rows_per_thread;
                    if (t < remaining_rows) {
                        rows_for_this_thread++;
                    }
                    
                    thread_data[t].end_row = std::min(current_row + rows_for_this_thread, n);
                    current_row = thread_data[t].end_row;
                    
                    global_thread_pool->enqueue([&thread_data, t]() {
                        find_pivot_row_partial_task(&thread_data[t]);
                    });
                }
                
                global_thread_pool->wait_all();
                
                for (int t = 0; t < num_threads; t++) {
                    if (thread_data[t].local_max_val > maxVal) {
                        maxVal = thread_data[t].local_max_val;
                        pivotRow = thread_data[t].local_pivot_row;
                    }
                }
            } else {
                for (int k = i + 1; k < n; k++) {
                    double val = fabs(A_copy[k * n + i]);
                    if (val > maxVal) {
                        maxVal = val;
                        pivotRow = k;
                    }
                }
            }
        }
        
        if (maxVal <= epsilon_d) {
            delete[] A_copy;
            delete[] b_copy;
            return false;
        }
        
        if (pivotRow != i) {
            for (int j = 0; j < n; j++) {
                std::swap(A_copy[i * n + j], A_copy[pivotRow * n + j]);
            }
            std::swap(b_copy[i], b_copy[pivotRow]);
        }
        
        double pivot = A_copy[i * n + i];
        for (int j = i; j < n; j++) {
            A_copy[i * n + j] /= pivot;
        }
        b_copy[i] /= pivot;
        
        if (n - i - 1 > 0) {
            if (use_multithreading && (n - i - 1 > 100)) {
                std::vector<ThreadData> thread_data(num_threads);
                
                int rows_per_thread = (n - i - 1) / num_threads;
                int remaining_rows = (n - i - 1) % num_threads;
                
                int current_row = i + 1;
                for (int t = 0; t < num_threads; t++) {
                    thread_data[t].A = A_copy;
                    thread_data[t].b = b_copy;
                    thread_data[t].n = n;
                    thread_data[t].current_pivot = i;
                    thread_data[t].start_row = current_row;
                    
                    int rows_for_this_thread = rows_per_thread;
                    if (t < remaining_rows) {
                        rows_for_this_thread++;
                    }
                    
                    thread_data[t].end_row = std::min(current_row + rows_for_this_thread, n);
                    current_row = thread_data[t].end_row;
                    
                    global_thread_pool->enqueue([&thread_data, t]() {
                        null_rows_task(&thread_data[t]);
                    });
                }
                
                global_thread_pool->wait_all();
            } else {
                for (int k = i + 1; k < n; k++) {
                    double factor = A_copy[k * n + i];
                    for (int j = i; j < n; j++) {
                        A_copy[k * n + j] -= A_copy[i * n + j] * factor;
                    }
                    b_copy[k] -= b_copy[i] * factor;
                }
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        
        if (n - i - 1 > 0) {
            if (use_multithreading && (n - i - 1 > 100)) {
                std::vector<BackSubThreadData> thread_data(num_threads);
                std::vector<double> partial_sums(num_threads, 0.0);
                
                int cols_per_thread = (n - i - 1) / num_threads;
                int remaining_cols = (n - i - 1) % num_threads;
                
                int current_col = i + 1;
                for (int t = 0; t < num_threads; t++) {
                    thread_data[t].A = A_copy;
                    thread_data[t].x = x;
                    thread_data[t].n = n;
                    thread_data[t].i = i;
                    thread_data[t].start_j = current_col;
                    
                    int cols_for_this_thread = cols_per_thread;
                    if (t < remaining_cols) {
                        cols_for_this_thread++;
                    }
                    
                    thread_data[t].end_j = std::min(current_col + cols_for_this_thread, n);
                    thread_data[t].partial_sum = &partial_sums[t];
                    current_col = thread_data[t].end_j;

                    global_thread_pool->enqueue([&thread_data, t]() {
                        compute_partial_sum_task(&thread_data[t]);
                    });
                }

                global_thread_pool->wait_all();
                
                for (int t = 0; t < num_threads; t++) {
                    sum += partial_sums[t];
                }
            } else {
                for (int j = i + 1; j < n; j++) {
                    sum += A_copy[i * n + j] * x[j];
                }
            }
        }
        
        x[i] = b_copy[i] - sum;
    }
    
    delete[] A_copy;
    delete[] b_copy;
    return true;
}

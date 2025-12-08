/* #ifndef SOLUTION_H
#define SOLUTION_H

#include <pthread.h>

bool solve(double* A, double* b, int n, double* x, int num_threads);

#endif */

#ifndef SOLUTION_H
#define SOLUTION_H

bool solve(double* A, double* b, int n, double* x, int num_threads);
void init_thread_pool(int num_threads);
void cleanup_thread_pool();

#endif


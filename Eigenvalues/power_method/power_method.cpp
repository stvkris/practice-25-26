#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "power_method.h"

void multiply_matrix_vector(double* A, double* x, double* result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += A[i * n + j] * x[j];
        }
    }
}

double norm(double* x, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return std::sqrt(sum);
}

void normalize(double* x, int n) {
    double nrm = norm(x, n);
    if (nrm == 0.0) nrm = 1.0;
    for (int i = 0; i < n; i++) {
        x[i] /= nrm;
    }
}

double product(double* x, double* y, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}

int power_method(double** A, int n, double* eigenvector,
                 double e, double* eigenvalue)
{
    double* x = new double[n];
    double* y = new double[n];
    
    for (int i = 0; i < n; i++) {
        x[i] = 1.0;
    }
    normalize(x, n);
    
    double lambda_old = 0.0;
    double lambda_new = 0.0;
    double error = e + 1.0;
    int max_iter = 1000;
    int iter = 0;
    
    for (iter = 0; iter < max_iter && error > e; iter++) {
        multiply_matrix_vector(*A, x, y, n);
        
        lambda_new = product(x, y, n);

        if (iter > 0) {
            error = std::abs(lambda_new - lambda_old);
        }

        for (int i = 0; i < n; i++) {
            x[i] = y[i];
        }
        normalize(x, n);
        
        lambda_old = lambda_new;
    }
    
    *eigenvalue = lambda_new;
    for (int i = 0; i < n; i++) {
        eigenvector[i] = x[i];
    }
    
    delete[] x;
    delete[] y;
    
    return 0;
}


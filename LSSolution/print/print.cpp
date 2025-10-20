#include <cstdlib>
#include <iostream>

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

double norm(double* A, double* x, double* b, int n)
{
    double epsilon_d = std::numeric_limits<double>::epsilon();
    
    double norm_top = 0.0;
    double norm_bot = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        norm_top += (sum - b[i]) * (sum - b[i]);
        norm_bot += b[i] * b[i];
    }
    if (norm_bot <= epsilon_d) {
        return -1.0;
    }
    return sqrt(norm_top / norm_bot);
}

#include <cmath>
#include <algorithm>

bool solve(double* A, double* b, int n, double* x)
{
    double epsilon_d = std::numeric_limits<double>::epsilon();
    
    for (int i = 0; i < n; i++) {
        
        int pivotRow = i;
        double maxVal = fabs(A[i * n + i]);
        for (int k = i + 1; k < n; k++) {
            double val = fabs(A[k * n + i]);
            if (val > maxVal) {
                maxVal = val;
                pivotRow = k;
            }
        }

        
        if (maxVal <= epsilon_d)
            return false;

        
        if (pivotRow != i) {
            for (int j = 0; j < n; j++)
                std::swap(A[i * n + j], A[pivotRow * n + j]);
            std::swap(b[i], b[pivotRow]);
        }

       
        double pivot = A[i * n + i];
        for (int j = i; j < n; j++)
            A[i * n + j] /= pivot;
        b[i] /= pivot;

        
        for (int k = i + 1; k < n; k++) {
            double factor = A[k * n + i];
            for (int j = i; j < n; j++)
                A[k * n + j] -= A[i * n + j] * factor;
            b[k] -= b[i] * factor;
        }
    }

    
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        x[i] = b[i] - sum;
    }

    return true;
}

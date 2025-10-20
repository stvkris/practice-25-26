#include <iostream>
#include <fstream>
#include <cmath>

bool from_file(char* filename, int n, double* A)
{
    std::ifstream file(filename);
    
    if (!file){
        return false;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j ++) {
            if (!(file >> A[i * n + j])) {
                return false;
            }
        }
    }
    
    return true;
}

bool by_formula(int k, int n, double* A)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j ++) {
            int ni = i+1;
            int nj = j+1;
            switch(k) {
                case 1:
                    A[i * n + j] = n - std::max(ni, nj) + 1;
                    break;
                case 2:
                    A[i * n + j] = std::max(ni, nj);
                    break;
                case 3:
                    A[i * n + j] = std::abs(ni - nj);
                    break;
                case 4:
                    A[i * n + j] = 1 / (ni + nj - 1);
                    break;
                default:
                    return false;
            }
        }
    }
    return true;
}

void init_b(double* b, double* A, int n)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < std::min(n, (int)std::lround((n + 1) / 2.0) + 1); j++) {
            sum += A[n * i + j];
        }
        b[i] = sum;
    }
}


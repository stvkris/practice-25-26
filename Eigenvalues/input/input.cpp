#include <iostream>
#include <fstream>
#include <cmath>

#include "input.h"

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
                    A[i * n + j] = n - std::max(i, j);
                    break;
                case 2:
                    if (i == j) {
                        A[i * n + j] = 2;
                    }
                    else if (abs(i - j) == 1) {
                        A[i * n + j] = 1;
                    }
                    else {
                        A[i * n + j] = 0;
                    }
                    break;
                case 3:
                    if (i == j && j + 1 < n) {
                        A[i * n + j] = 1;
                    }
                    else if (j + 1 == n) {
                        A[i * n + j] = i + 1;
                    }
                    else if (i + 1 == n) {
                        A[i * n + j] = j + 1;
                    }
                    else {
                        A[i * n + j] = 0;
                    }
                    break;
                case 4:
                    A[i * n + j] = (double)(1) / (i + j + 1);
                    break;
                default:
                    return false;
            }
        }
    }
    return true;
}

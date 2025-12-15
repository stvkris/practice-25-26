#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "print.h"

void print_m(double* A, int n, int m)
{
    int s = std::min(n, m);
    std::cout << std::scientific << std::setprecision(3);
        
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < s; j++) {
                std::cout << " " << std::setw(10) << A[i * n + j];
            }
            std::cout << std::endl;
        }
}

void print_v(double* b, int n, int m)
{
    int s = std::min(n, m);
    for (int i = 0; i < s; i++) {
        std::cout << b[i] << "\n";
    }
}

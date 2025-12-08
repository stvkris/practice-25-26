#include <iostream>
#include <cstdlib>
#include <chrono>
#include "input.h"
#include "print.h"
#include "solution.h"

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Incorrect characteristics";
        return 1;
    }
    
    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);
    int num_threads = std::atoi(argv[4]);
    
    if (num_threads <= 0) {
        std::cerr << "Number of threads must be positive\n";
        return 1;
    }
    
    double* A = new double[n * n];
    
    bool flag;
    
    if (k < 5 && k > -1) {
        if (k == 0) {
            if (argc != 5) {
                std::cerr << "Where is filename?!";
                delete[] A;
                return 1;
            }
            char* filename = argv[5];
            flag = from_file(filename, n, A);
        } else {
            flag = by_formula(k, n, A);
        }
    } else {
        std::cerr << "Incorrect value of k";
        delete[] A;
        return 1;
    }
    
    if (!flag) {
        std::cerr << "Error with initialization";
        delete[] A;
        return 1;
    }
    
    double* b = new double [n];
    double* x = new double [n];
    
    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
    }
    
    init_b(b, A, n);
    
    std::cout << "Matrix A:\n";
    print_m(A, n, m);
    std::cout << "Vector b:\n";
    print_v(b, n, m);
    
    double* A_original = new double[n * n];
    double* b_original = new double[n];

    for (int i = 0; i < n * n; i++) {
        A_original[i] = A[i];
    }
    for (int i = 0; i < n; i++) {
        b_original[i] = b[i];
    }
    
    if (num_threads > n) {
        num_threads = n;
    }
    
    auto start1 = std::chrono::high_resolution_clock::now();
    
    flag = solve(A, b, n, x, num_threads);
    
    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
    
    if (!flag) {
        std::cerr << "Error with finding solution\n";
    }
    else {
        std::cout << "Solution x:\n";
        print_v(x, n, m);
        
        auto start2 = std::chrono::high_resolution_clock::now();
        
        double sol_norm = norm(A_original, x, b_original, n, num_threads);
        
        auto end2 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
        
        std::cout << "Norm = " << sol_norm << "\n";
        
        double total_seconds = (duration1.count() + duration2.count()) / 1000.0;
        std::cout << "Time: " << total_seconds << " sec\n";
    }
    
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] A_original;
    delete[] b_original;
    
    return 0;
}


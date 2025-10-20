#include <iostream>
#include <cstdlib>
#include <ctime>
#include "input.h"
#include "print.h"
#include "solution.h"

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Incorrect characteristics";
        return 1;
    }
    
    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);
    
    double* A = new double[n * n];
    
    bool flag;
    
    if (k < 5 && k > -1) {
        
        if (k == 0) {
            
            if (argc != 5) {
                std::cerr << "Where is filename?!";
                delete[] A;
                return 1;
            }
            char* filename = argv[4];
            
            flag = from_file(filename, n, A);
            
            if (!flag) {
                std::cerr << "Error with reading file";
                delete[] A;
                return 1;
            }
        }
        
        else {
            
            flag = by_formula(k, n, A);
            
            if (!flag) {
                std::cerr << "Error with inizialition by formula";
                delete[] A;
                return 1;
            }
        }
    }
    
    else {
        std::cerr << "Incorrect value of k";
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
    std::cout << "Vector b;\n";
    print_v(b, n, m);
    
    double* A_original = new double[n * n];
    double* b_original = new double[n];

    for (int i = 0; i < n * n; i++) {
        A_original[i] = A[i];
    }
    for (int i = 0; i < n; i++) {
        b_original[i] = b[i];
    }
    
    clock_t start = clock();
    
    flag = solve(A, b, n, x);
    
    clock_t end = clock();
    
    if (!flag) {
        std::cerr << "Error with finding solution\n";
    }
    else {
        std::cout << "Solution x:\n";
        print_v(x, n, m);
        
        double sol_norm = norm(A_original, x, b_original, n);
        std::cout << "Norm = " << sol_norm << "\n";
    }
    
    std::cout << "Time: " << double(end - start) / CLOCKS_PER_SEC << " sec/n \n";
    
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] A_original;
    delete[] b_original;
    
    return 0;
}

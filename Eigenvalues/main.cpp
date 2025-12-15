#include <iostream>
#include <cstdlib>
#include <ctime>
#include "input.h"
#include "print.h"
#include "power_method.h"

int main(int argc, char * argv[]) {
    
    if (argc < 5 || argc > 6) {
        std::cerr << "Incorrect characteristics";
        return 1;
    }
    
    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    double e = std::atof(argv[3]);
    int k = std::atoi(argv[4]);
    
    double* A = new double[n * n];
    
    bool flag;
    
    if (k < 6 && k > -1) {
        
        if (k == 0) {
            
            if (argc != 6) {
                std::cerr << "Where is filename?!";
                delete[] A;
                return 1;
            }
            char* filename = argv[5];
            
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
    
    std::cout << "Matrix A:\n";
    print_m(A, n, m);
    
    double* eigenvector = new double[n];
    double eigenvalue;
    
    clock_t start = clock();
    power_method(&A, n, eigenvector, e, &eigenvalue);
    clock_t end = clock();
    
    std::cout << "Eigenvector:\n";
    print_v(eigenvector, n, m);
    
    std::cout << "Maximun of eigenvalues: " << eigenvalue << "\n";
    
    std::cout << "Time: " << double(end - start) / CLOCKS_PER_SEC << " sec \n";
    return 0;
}

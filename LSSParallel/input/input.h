#ifndef INPUT_H
#define INPUT_H

bool from_file(char* filename, int n, double* A);
void init_b(double* b, double* A, int n);
bool by_formula(int k, int n, double* A);

#endif

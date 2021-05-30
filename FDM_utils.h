#ifndef FDM_UTILS_H
#define FDM_UTILS_H

double* thomas_method(int n, double *a, double *b);
double* sor_method(int n, double *a, double *b, double relax, int max_iter);

#endif 
#ifndef HEADER_H
#define HEADER_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


double f(double x);
double* single_division_scheme(double **A, int N);
double Legandr(int n, double x);
double func(int N, double x, int k1, int k2, double a, double b);
double trapeze(int k1, int k2, double a, double b, int n, int N);
double Recalculation_Integral(int k1, int k2, double a, double b, int N);
double * Ak(double a, double b, int N);
double Polinom(double x, double a, double b, int n, double *A);
int N_Polinom(double eps, double a, double b);

#endif //header
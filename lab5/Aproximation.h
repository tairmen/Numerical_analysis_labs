#pragma once
//created by Tair Abduraimov
//kv-52
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>


double func(double x);
double Legendre(double x, int m);
double trapezium(double a, double b, int n, int i, int k);
double* single_division(double **A, int n);
double** fill_matrix(int n, double a, double b);
double aprox(int n, double x);
int eps_step(double eps, double a, double b);

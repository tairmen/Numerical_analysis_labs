#include "header.h"


double f(double x){
	double p = (double)1 / 3;
	return 0.2*pow(M_E, pow(x, p))*log(x)*sin(3*x);
}

double* single_division_scheme(double **A, int N) {
	int i, j, k, l;
	double tmp1, tmp2, tmp3;
	double* res = new double[N];

	for (j = 0; j < N; j++) {
		tmp1 = A[j][j];
		for (i = j; i < N + 1; i++) {
			A[j][i] /= tmp1;
		}
		for (k = j + 1; k < N; k++) {
			tmp2 = A[k][j];
			for (l = j; l < N + 1; l++) {
				A[k][l] -= A[j][l] * tmp2;
			}
		}
	}
	for (j = N - 1; j >= 0; j--) {
		tmp3 = 0;
		for (i = N - 1; i > j; i--) {
			tmp3 += A[j][i] * res[i];
		}
		res[j] = A[j][N] - tmp3;
	}
	return res;
}

double Legandr(int n, double x){
	double Ln1, Ln = x, Ln_1 = 1;
	if (n == 0) return Ln_1;
	if (n == 1) return Ln;

	int i = 1;
	while (i < n){
		Ln1 = (1.0 * (2 * i + 1) / (i + 1)) * x * Ln - (1.0 * i / (i + 1)) * Ln_1;
		Ln_1 = Ln;
		Ln = Ln1;
		++i;
	}
	return Ln;
}

double func(int N, double x, int k1, int k2, double a, double b){
	double t = (2 * x - a - b) / (b - a);
	if (k2 == N)
		return  f(x)*Legandr(k1, t);
	return  Legandr(k1, t)*Legandr(k2, t);
}

double trapeze(int k1, int k2, double a, double b, int n, int N){
	double h, x, integral = 0;
	h = (b - a) / n;
	for (x = a; x < b; x += h){
		integral += h*(func(N, x, k1, k2, a, b) + func(N, x + h, k1, k2, a, b)) / 2;
	}
	return integral;
}

double Recalculation_Integral(int k1, int k2, double a, double b, int N){
	double error, In, I2n, eps=1e-6;
	int n = (int)ceil((b - a) / sqrt(sqrt(eps)));
	In = trapeze(k1, k2, a, b, n, N);
	I2n = trapeze(k1, k2, a, b, 2 * n, N);
	error = fabs(In - I2n) / 15;

	while (error > eps) {
		In = I2n;
		n *= 2;
		I2n = trapeze(k1, k2, a, b, n, N);
		error = fabs(In - I2n) / 15;
	}
	return I2n;
}

double * Ak(double a, double b, int N){
	double *result = new double[N];
	double ** A = new double*[N];
	for (int i = 0; i < N + 1; i++)
		A[i] = new double[N + 1];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N + 1; j++) {
			A[i][j] = A[j][i] = Recalculation_Integral(i, j, a, b, N);
		}
		A[i][N] = Recalculation_Integral(i, N, a, b, N);
	}

	result = single_division_scheme(A, N);
	return result;
}

double Polinom(double x, double a, double b, int n, double *A){
	double  t = (2 * x - b - a) / (b - a);
	double function = 0;
	for (int i = 0; i < n; i++)
		function += A[i] * Legandr(i, t);
	return function;
}

int N_Polinom(double eps, double a, double b){
	double y, x, Error, step = 0.2;
	int n = 2, k = (int)(b - a) / step;

	do{
		n *= 2;
		Error = 0;
		x = a;
		for (int i = 0; i <= k; ++i){
			double *A = Ak(a, b, n);
			y = f(x) - Polinom(x, a, b, n, A);
			Error += y*y;
			x += step;
		}
		Error = sqrt(Error / (k + 1));
	} while (Error > eps);

	return n;
}


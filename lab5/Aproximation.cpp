#include "Aproximation.h"

double func(double x) {
	return (0.1*x*x*exp(pow(x,1.0/3))*sin(3*x));
}

double Legendre(double x, int m) {
	double Pn, P, tmp;
	int i;
	i = 2;
	P = 1;
	Pn = x;
	if (m == 0) {
		return P;
	}
	else if (m == 1) {
		return Pn;
	}
	else {
		while (i <= m) {
			tmp = Pn;
			Pn = x*Pn*(2 * i + 1) / (i + 1) - i*P / (i + 1);
			P = tmp;
			i++;
		}
		return Pn;
	}
}

double trapezium(double a, double b, int n, int i, int k) {
	double h, x;
	int p;
	double integral = 0;
	h = (b - a) / n;
	x = a;
	if (i < 0) {
		for (p = 0; p < n; p++) {
			integral += h*(func(x) * Legendre(x, k) + func(x + h) * Legendre(x + h, k)) / 2;
			x += h;
		}
	}
	else {
		for (p = 0; p < n; p++) {
			integral += h*(Legendre(x, i) * Legendre(x, k) + Legendre(x + h, i) * Legendre(x + h, k)) / 2;
			x += h;
		}
	}
	return integral;
}

double Runge(double a, double b, double eps, int n, int i, int k) {
	
		double Rn, In, I2n;
		In = trapezium(a, b, n, i, k);
		n = n * 2;
		I2n = trapezium(a, b, n, i, k);
		Rn = fabs((In - I2n) / In);
		while (Rn > eps) {
			In = I2n;
			n = n * 2;
			I2n = trapezium(a, b, n, i, k);
			Rn = fabs((In - I2n) / In);
		}
		return I2n;
	
}

double** fill_matrix(int n, double a, double b) {
	int i, k;
	double** matrix;
	matrix = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		matrix[i] = (double *)malloc((n + 1) * sizeof(double));
	for (k = 0; k < n; k++) {
		for (i = 0; i < n; i++) {
			matrix[k][i] = Runge(a, b, 0.01, 10, i, k);//trapezium(a, b, 200, i, k);
		}
		matrix[k][n] = Runge(a, b, 0.01, 10, -1, k);//trapezium(a, b, 200, -1, k);
	}
	return matrix;
}

double* single_division(double **A, int n) {
	double* res;
	res = (double *)malloc(n * sizeof(double));
	int i, j, k, l;
	double tmp1, tmp2, tmp3;
	for (j = 0; j < n; j++) {
		tmp1 = A[j][j];
		for (i = j; i < n+ 1; i++) {
			A[j][i] /= tmp1;
		}
		for (k = j + 1; k < n; k++) {
			tmp2 = A[k][j];
			for (l = j; l < n + 1; l++) {
				A[k][l] -= A[j][l] * tmp2;
			}
		}
	}
	for (j = n - 1; j >= 0; j--) {
		tmp3 = 0;
		for (i = n - 1; i > j; i--) {
			tmp3 += A[j][i] * res[i];
		}
		res[j] = A[j][n] - tmp3;
	}
	return res;
}

double* Gauss(double **m, int n) {
	int i, j, h;
	double tmp;
	double* res;
	double** matrix;
	matrix = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		matrix[i] = (double *)malloc((n + 1) * sizeof(double));
	for (i = 0; i < n; i++)
		for (j = 0; j < n + 1; j++)
			matrix[i][j] = m[i][j];
	res = (double *)malloc(n * sizeof(double));
	for (i = 0; i < n; i++) {
		tmp = matrix[i][i];
		for (j = 0; j < n + 1; j++) {
			//print_matrix(matrix, n);		
			matrix[i][j] = matrix[i][j] / tmp;
		}
		for (h = 0; h < n; h++) {
			tmp = matrix[h][i];
			//print_matrix(matrix, n);
			for (j = 0; j < n + 1; j++) {
				if (h == i)
					break;
				matrix[h][j] = matrix[h][j] - matrix[i][j] * tmp;
			}
		}
	}
	for (i = 0; i < n; i++) {
		res[i] = matrix[i][n];
	}
	return res;
}

double aprox(int n, double x) {
	int i, j, a, b;
	double Apr = 0;
	double* A;
	double** matrix;
	A = (double *)malloc(n * sizeof(double));
	matrix = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		matrix[i] = (double *)malloc((n + 1) * sizeof(double));
	//pp = eps_step(0.1, 0, 1);
	//printf("%d ", pp);
	a = (int)x;
	b = a + 1;
	matrix = fill_matrix(n, a, b);
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n + 1; j++) {
	//		printf("%20.12f\n", matrix[i][j]);
	//	}
	//}
	//printf("\n");
	A = single_division(matrix, n);
	//for (i = 0; i < n; i++) {
	//	printf("%20.12f\n", A[i]);
	//}
	for (i = 0; i < n; i++) {
		Apr += A[i] * Legendre(x, i);
	}
	return Apr;
}

int eps_step(double eps, double a, double b) {
	int i = 0;
	double h, x;
	double integral;
	int n = 40;
	h = (b - a) / n;
	do {
		integral = 0;
		for (x = a; x <= b; x += h) {
			integral += h*((func((x + h) / 2) - aprox(i, (x+h)/2)) * (func((x + h) / 2) - aprox(i, (x + h) / 2)));
		}
		i++;
		//printf("int = %e\n", integral);
		//printf("pox = %e\n", func(x) - aprox(i, x));
	} while (sqrt(fabs(integral) / (b - a)) > eps);
	return i;
}
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

double a = 20;
double b = 30;

double f(double x){
	return pow(M_E, pow((x - 15), 1.0 / 4))*sin(1.2*x);
}

double **Make_SLAR(int N){
	double **res = new double*[N];
	for (int i = 0; i < N; i++)
		res[i] = new double[N + 1];

	double h = (b - a) / N;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N + 1; j++)
			res[i][j] = 0;

	for (int i = 0; i < N; i++){
		res[i][i] = 4 * h;
		if (i>0)
			res[i][i - 1] = h;
		if (i < N-1)
			res[i][i + 1] = h;
		double x_prev = a + h * (i - 1);
		double x_i = a + h * i;
		double x_next = a + h * (i + 1);
		res[i][N] = 6 * (f(x_next) - 2 * f(x_i) + f(x_prev)) / h;
	}
	return res;
}


double * Gaus_Jordan(int N){
	double **A = Make_SLAR(N);
	double *res = new double[N];
	double coeficient;

	for (int k = 0; k < N; k++){

		coeficient = A[k][k];
		for (int j = 0; j < N + 1; j++)
			A[k][j] = A[k][j] / coeficient;

		for (int i = 0; i < N; i++)
			if (i != k){
				coeficient = A[i][k];
				for (int j = 0; j < N + 1; j++)
					A[i][j] -= A[k][j] * coeficient;
		}
	}

	for (int i = 0; i < N; i++)
		res[i] = A[i][N];

	return res;
}


double* findAi(int N){
	double* A = new double[N];
	double h = (b - a) / N;
	double x = a;
	for (int i = 0; i < N; ++i)	{
		A[i] = f(x);
		x += h;
	}
	return A;
}

double* findDi(double* C, int N){
	double* D = new double[N];
	double h = (b - a) / N;

	D[0] = C[0];
	for (int i = 1; i < N; ++i)
		D[i] = (C[i] - C[i - 1]) / h;
	return D;
}

double* findBi(double* C, double* D, int N){
	double* B = new double[N];
	double h = (b - a) / N;
	B[0] = 0;

	for (int i = 1; i < N; ++i)
		B[i] = (h * C[i] / 2) - (h * h * D[i] / 2) + (f(a+h*i) - f(a+h*(i-1))) / h;
	return B;
}

double Spline(double x, double* A, double* B, double* C, double* D, int N){
	double h = (b - a) / N;
	int i = 0;

	//Номер відрізку
	while (x >= (a + i * h))
		i++; 
	i--;

	double xi = a + i * h;
	//Значення сплайну в точці
	return A[i] + B[i] * (x - xi) + C[i] * (x - xi) * (x - xi) / 2 + D[i] * (x - xi) * (x - xi) * (x - xi) / 6;
}

void SplineInterpolation(int N){

	double* A;
	double* C;
	double* D;
	double* B;

	A = findAi(N);
	C = Gaus_Jordan(N);
	D = findDi(C, N);
	B = findBi(C, D, N);

	ofstream fout("results.txt");
	int k = (b - a) * 20+1; //кількість точок, які обчислюються на заданому інтервалі
	double step = (b - a) / (double)k;
	double x = a, y;

	for (int i = 0; i <= k; i++){
		y = Spline(x, A, B, C, D, N);
		fout << x;
		fout << " ; "; 
		fout << y;
		fout << endl;
		x += step;
	}
	fout.close();
	return;
}


int main(){
	int N = 50;
	printf("Building spline with N = %d\n", N);
	SplineInterpolation(N);
	printf("Building finished \n");

	return 0;
}
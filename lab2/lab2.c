#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//double P = M_PI;

double  myfunction(double);
double  frst_derivative(double);
int     roots_selection(double, double, double, double*);
int*  iteration(double*, double, int);
int* bisection(double*, double, int);

main() {
	double mas[5];
	int *n_it, *n_b;
	int i, k;
	double step = 0.2;//step < 0.5
	double p = 1e+1;

	k = roots_selection(-10, 10, step, mas);
	printf("%d\n", k);

	n_it = iteration(mas, step, k);
	printf("------------------------------------------------------\n\n");
	n_b = bisection(mas, step, k);


	printf("------------------------------------------------------\n\n");
	printf(" EPS     ITERATION    BISECTION \n\n");
	for (i = 0; i <= 4; i++)
		printf("%0.E     %2.0d           %d\n", p /= 1000, n_it[i], n_b[i]);

	system("PAUSE");

	return;
}

double myfunction(double x) {
	return sin(x*x) - x/ 2 + 0.5;
}

double frst_derivative(double x) {
	return 2 * x*cos(x*x) - 0.5;
}

int roots_selection(double border1, double border2, double step, double *mas) {

	int i = -1;
	double a;

	for (a = border1; a <= border2; a += step)
		if (myfunction(a) * myfunction(a + step) < 0)
			*(mas + (++i)) = a;
	return i;
}

int* iteration(double *mas, double step, int sum_roots) {
	double m, M, a, b;
	double alfa, q, eps;
	double x1, x;
	int n, i, t = 0;

	int* n_it = (int*)malloc(sizeof(int) * 4);

	for (i = 0; i <= sum_roots; i++){

		printf(" EPS            ROOT_%d              ACCURATELY\n\n", i + 1);

		a = *(mas + i);
		b = a + step;


		if (fabs(frst_derivative(a)) > fabs(frst_derivative(b))) {
			M = frst_derivative(a);
			m = frst_derivative(b);
		}
		else {
			M = frst_derivative(b);
			m = frst_derivative(a);
		}
		alfa = 1 / M;
		q = 1 - fabs(m / M);

		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
			x = (b + a) / 2;
			x1 = x;
			n = 0;
			do {
				n++;
				x = x1;
				x1 = x - alfa*myfunction(x);
			} while (fabs(x1 - x) > (1 - q) / q*eps);
			if (i == 0)
				*(n_it + t++) = n;
			printf("%0.E  %20.15f  %20.15f   \n", eps, x1, fabs(fabs(x1) - fabs(x))*q / (1 - q));

		}
		printf("\n");

	}
	return n_it;
}

int* bisection(double *mas, double step, int sum_roots){
	double a, b, c;
	double x, eps;
	int i, n, t=0;

	int* n_b = (int*)malloc(sizeof(int) * 4);

	for (i = 0; i <= sum_roots; i++){
		printf(" EPS            ROOT_%d              ACCURATELY\n\n", i + 1);

		a = mas[i];
		b = a + step;

		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3){
			n = 0;
			do{
				c = (a + b) / 2;
				if (myfunction(a)*myfunction(c) < 0)
					b = c;
				if (myfunction(b)*myfunction(c) < 0)
					a = c;
				n += 1;
			} while (fabs(b - a) >= 2 * eps);
			x = (a + b) / 2;
			if (i == 0)
				*(n_b + t++) = n;
			printf("%0.E  %20.15f  %20.15f   \n", eps, x, fabs((b-a)/2));
		}
	}
	return n_b;
}


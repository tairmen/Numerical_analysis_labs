#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES 



double my_func(double x){
	return 1 / x;
}

double primitive_func(double x){
	return log(x);
}

double integrate(double a, double b, int n){
	double h, x, integral = 0;
	h = (b - a) / n;
	for (x = a; x < b; x += h){
		integral += h*(my_func(x) + my_func(x + h)) / 2;
	}
	return integral;
}

double runge(double a, double b, double eps){
	double r, In, I2n;
	int n = (int)ceil((b-a)/ sqrt(eps));

	In = integrate(a, b, n);
	I2n = integrate(a, b, 2 * n);
	r = fabs(In - I2n) / 3;

	while (r > eps) {
		In = I2n;
		n *= 2;
		I2n = integrate(a, b, n);
		r = fabs(In - I2n) / 3;
	}

	return I2n;
}


int main(){
	double a, b, h, error, value, integral, eps = 1e-2;
	double n;
	a = 1;
	b = 375;
	h = 0.013;//h визначено аналітично
	n = (int)ceil((b - a) / h);
	value = primitive_func(b) - primitive_func(a);
	error = fabs(integrate(a, b, n) - value);

	printf("Task1: \n");
	printf("Epsilon     Step    Value    Error\n");
	printf("%8.6f  %8.6f  %5.3f  %10.8f\n", eps, (b-a) / n, value, error);
	printf("\n\n");


	n = (int)ceil((b-a)/ sqrt(error));
	h = (b - a) / n;
	printf("Task2: \n");
	integral = runge(a, b, error);
	printf("Integral = %f \n", integral);
	printf("Epsilon       Step      Error\n");
	printf("%10.8f  %8.6f  %10.8f\n", error, h, fabs(integral - integrate(a, b, n)));
	return;
}
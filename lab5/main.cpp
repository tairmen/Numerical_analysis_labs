#include "Aproximation.h"

int main() {
	//printf("%e", func(2));
	int n;
	int a = 0;
	int b = 6;
	double x , res1, res2;
	n = eps_step(0.1, a, b);
	printf("n=%d\n", n);
	FILE *file;
	file = fopen("table.csv", "wb");
	fprintf(file, "%4s;%24s;%24s\n", "x", "aproximation", "native function");
	for (x = a; x <= b; x += 0.1) {
		res1 = aprox(n, x);
		res2 = func(x);
		fprintf(file, "%4.2f;%24.17f;%24.17f\n", x, res1, res2);
	}
	fclose(file);
	return 0;
}
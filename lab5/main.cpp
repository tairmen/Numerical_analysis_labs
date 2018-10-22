#include "header.h"

int main(){
	ofstream  stream;
	double a = 2, b = 9;
	double eps = 0.1;
	cout << "calculation..." << endl;

	int N = N_Polinom(eps, a, b);
	cout << "N = " << N << endl;
	stream.open("result.txt");
	for (double x = a; x <= b; x += 0.2){
		double *A = Ak(a, b, N);
		double res = Polinom(x, a, b, N, A);
		stream << setw(10) << fixed << x << ";"<< res << endl;
		printf("%5.3f  %10.8f\n", x, res);
		
	}
	stream.close();

	return 0;
}
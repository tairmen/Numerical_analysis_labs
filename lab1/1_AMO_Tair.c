//copyright by Tair Abduraimov
//group KV-52; project - calculate function ln(x)
#include <stdio.h>
#include <math.h>

#define M_LN2          0.69314718055994530942

void makloren(double eps, double x, int p){

	double Lk = 1.0;
	double L;
	int m = 1;
	double xx;
	double z;
	double a;
	int k = 0;
	double R, abs_delta;
	
	xx = x;
	while (xx >= 2){
		xx = xx/2;
		m = m + 1;
	}	
	//printf("%15d\n",m);
	z = x/(pow(2,m));
	//printf("%15e\n",z);
	a = (1-z)/(1+z);
	L = m*M_LN2;
	//printf("%15e\n",Lk);
	ch = x;
	sum = x;
	while (ch > eps){
		//Lk = pow(a,2*k-1) / (2*k-1);
		//ch = (-1) ^ k*x ^ (2k + 1) / (2k + 1)!;

		//printf("%15e\n",Lk);
		//L -= 2*Lk;
		k = k + 2;
		ch = -ch*x ^ 2 / k*(k + 1);
		sum += ch;
	}
	R = ;

	abs_delta = fabs(L - sin(x));

	if (p==0)
		printf("|%20e |%6d |%20e |%20e |\n", eps, k, abs_delta, R);
	if (p==1)
		printf("|%10f |%20e |%20e |\n", x, abs_delta, R);
}


int main(){
	double a = 0.98;
	double b = 5.5;
	double x =(a + b) / 2;
	double h = fabs((b - a) / 10);
	double eps;
	int i;
	printf("for x = %f  :\n", x);
	printf("---------------------------------------------------------------------------\n");
	printf("|%20s |%6s |%20s |%20s |\n", "eps", "n", "abs_pogr", "R");
	printf("---------------------------------------------------------------------------\n");

	eps = 1e-2;
	while (eps>1e-14){
		makloren(eps, x, 0);
		eps = eps/1000;
	}
	
	eps = 1e-8;
	printf("\n\nfor eps = %e  :\n",eps);
	printf("---------------------------------------------------------\n");
	printf("|%10s |%20s |%20s |\n", "X", "abs_pogr", "R");
	printf("---------------------------------------------------------\n");

	x = a;
	for (i = 1; i <= 11; i++){
		makloren(eps, x, 1);
		x = a + h*i;
	}

	return 0;
}
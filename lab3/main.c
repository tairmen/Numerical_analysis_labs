#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// #define N 3
//#define M 4

void print_matrix(double **matrix, int n, int m){
	int i, j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			printf("%8.2f", matrix [i][j]);
		}
		printf("\n");
	}
	return;
}

double* Task1(double **matrix, int n, int m){
	int i, j, line, row, k, line_new=0;
	double *M, *res, **new_matrix, R, max;
	M = (double *)malloc(n * sizeof(double));
	res = (double *)malloc(n * sizeof(double));
	new_matrix = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		new_matrix[i] = (double *)malloc(m * sizeof(double));

	for (k = 0; k < n - 1; k++){

		max = 0;
		for (i = 0; i < n; i++)
			for (j = 0; j < m - 1; j++)
				if (fabs(matrix[i][j]) > max){
					max = fabs(matrix[i][j]);
					line = i;
					row = j;
				}

		for (j = 0; j < m; j++)
			new_matrix[line_new][j] = matrix[line][j];
		line_new++;

		for (i = 0; i < n; i++)
			M[i] = -(matrix[i][row] / max);

		for (i = 0; i < n; i++)
			for (j = 0; j < m; j++)
				if (i != line)
					matrix[i][j] += matrix[line][j] * M[i];

		for (j = 0; j < m; j++)
			matrix[line][j] = 0;

		for (i = 0; i < n; i++)
			matrix[i][row] = 0;

	}
	print_matrix(matrix, n, m);
	for (i = 0; i < n; i++)
		if (matrix[i][m - 1] != 0)
			line = i;

	for (j = 0; j < m; j++)
		new_matrix[line_new][j] = matrix[line][j];

	for (j = 0; j < m-1; j++)
		res[j] = 0;

	for (i = n - 1; i >= 0; i--){
		R = new_matrix[i][m - 1];
		for (j = 0; j < m - 1; j++)
			if ((new_matrix[i][j] != 0) && (res[j] != 0))
				R -= new_matrix[i][j] * res[j];
		for (j = 0; j < m - 1; j++)
			if ((new_matrix[i][j] != 0) && (res[j] == 0))
				res[j] = R / new_matrix[i][j];
	}

	return res;
}

double* Task2(double **matrix, int n, int m, double eps){
	double *res, *help_res, max, q, m_norma_x;
	int i, j, k;
	res = (double *)malloc(n * sizeof(double));
	help_res = (double *)malloc(n * sizeof(double));

	for (i = 0; i < n; i++){
		max = matrix[i][i];
		for (j = 0; j < m; j++)
			matrix[i][j] = matrix[i][j] / max;
	}

	q = 0;
	for (i = 0; i < n; i++){
		max = 0;
		for (j = 0; j < m - 1; j++)
			if (i != j)
				max += fabs(matrix[i][j]);
		if (max > q)
			q = max;
	}

	for (i = 0; i < n; i++)
		help_res[i] = matrix[i][m-1];

	for (i = 0; i < n; i++){
		res[i] = matrix[i][m-1];
		for (j = 0; j < m - 1; j++)
			if (i != j)
				res[i] -= matrix[i][j] * help_res[j];
	}

	for (i = 0; i < n; i++)
		help_res[i] = res[i];

	k = 2;
	do{
		for (i = 0; i < n; i++){
			res[i] = matrix[i][m - 1];
			for (j = 0; j < m - 1; j++)
				if (i != j)
					res[i] -= matrix[i][j] * help_res[j];
		}

		for (i = 0; i < n; i++)
			help_res[i] = res[i] - help_res[i];

		m_norma_x = fabs(help_res[0]);
		for (i = 1; i < n; i++)
			if (fabs(help_res[i])>m_norma_x)
				m_norma_x = fabs(help_res[i]);

		for (i = 0; i < n; i++)
			help_res[i] = res[i];

		k++;
	} while (m_norma_x > eps*((1 - q) / q));

	printf("Number of iteration= %d\n", k);

	return res;
}


int main(){
	double **My_matrix, **Modern_matrix, *task1_matrix, *task2_matrix, eps;
	int i, m = 5, n = 4;

	My_matrix = (double **)malloc(n * sizeof(double *));
	for (i  = 0; i < n; i++)
		My_matrix[i] = (double *)malloc(m * sizeof(double));

	Modern_matrix = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		Modern_matrix[i] = (double *)malloc(m * sizeof(double));
	
	My_matrix[0][0] = 8; My_matrix[0][1] = 4; My_matrix[0][2] = 8; My_matrix[0][3] = 20; My_matrix[0][4] = 148;
	My_matrix[1][0] = 8; My_matrix[1][1] = 27; My_matrix[1][2] = 12; My_matrix[1][3] = 6; My_matrix[1][4] = 87;
	My_matrix[2][0] = 16; My_matrix[2][1] = 13; My_matrix[2][2] = 46; My_matrix[2][3] = 16; My_matrix[2][4] = 157;
	My_matrix[3][0] = 19; My_matrix[3][1] = 10; My_matrix[3][2] = 7; My_matrix[3][3] = 17; My_matrix[3][4] = 169;

	Modern_matrix[0][0] = 50; Modern_matrix[0][1] = -23; Modern_matrix[0][2] = 0; Modern_matrix[0][3] = -6; Modern_matrix[0][4] = 91;
	Modern_matrix[1][0] = 8; Modern_matrix[1][1] = 27; Modern_matrix[1][2] = 12; Modern_matrix[1][3] = 6; Modern_matrix[1][4] = 87;
	Modern_matrix[2][0] = 16; Modern_matrix[2][1] = 13; Modern_matrix[2][2] = 46; Modern_matrix[2][3] = 16; Modern_matrix[2][4] = 157;
	Modern_matrix[3][0] = 8; Modern_matrix[3][1] = -19; Modern_matrix[3][2] = 4; Modern_matrix[3][3] = 34; Modern_matrix[3][4] = 209;
	print_matrix(My_matrix, n, m);
	task1_matrix = Task1(My_matrix, n, m);
	printf("Roots Task1:\n\n");
	for (i = 0; i < m - 1; i++)
		printf("X(%d)=%f\n", i, task1_matrix[i]);
	
	printf("\n\n");
	eps = 1e-4;
	printf("Task2:\n\n");
	task2_matrix = Task2(Modern_matrix, n, m, eps);
	printf("EPS=%f\n\n", eps);
	printf("Roots Task2:\n\n");
	for (i = 0; i < m - 1; i++)
		printf("X(%d)=%f\n", i, task2_matrix[i]);
	
	return;
}
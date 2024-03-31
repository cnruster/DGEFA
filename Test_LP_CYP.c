// Test utilities for
// DGEFA+DGESL functions for Gaussian elimination linear equation solver
// Translated into C by Yiping Cheng, Beijing Jiaotong University, China
// Email: ypcheng@bjtu.edu.cn

#include "blas.h"
#include "blas_cyp.h"
#include <stdio.h>
#include <assert.h>
#include <time.h>

typedef unsigned int UINT;
void ILASEED(UINT sid);
UINT ILARAN();
void test_LP();
void test_CYP();
void compare_LP_CYP_times();

int main()
{
	// generate different seed with different time
	ILASEED(time(NULL));

	puts("% Testing correctness of my translation of DGEFA,DGESL,DGESLT of LinPack in C\n");
	test_LP();

	puts("% Testing correctness of CYP's implementation of DGEFA,DGESL,DGESLT in C\n");
	test_CYP();

	puts("% Please copy the above outputs and paste them into Matlab to verify correctness\n");

	puts("% Now comparing time performances of LinPack and CYP's implementations");
	puts("% It may take tens of minutes. Press Ctrl+C to exit");
	puts("% The results will show LinPack is better than CYP");
	compare_LP_CYP_times();
}



int rand_int(int min, int max)
{
	assert(min <= max);
	UINT rand_num = ILARAN(); // rand_num being unsigned is required
	return  min + (rand_num % (UINT)(max - min + 1));
}

double rand_double()
{
	return (double)rand_int(-5, 5) / (double)rand_int(1, 10);
}


void rand_vector(double* v, int n)
{
	for (int i = 0; i < n; i++) {
		v[i] = rand_double();
	}
}

void rand_matrix(double* A, int lda, int n)
{
	int i, j;

	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
			A[i+lda*j] = rand_double();
		}
	}
}

void print_vector(double* v, int n, char* name)
{
	printf("%s=[", name);
	for (int i = 0; i < n; i++) {
		printf("%.18f; ", v[i]);
	}
	printf("];\n");
}

void print_matrix(double* A, int lda, int n, char* name)
{
	int i, j;

	printf("%s=[", name);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%.18f, ", A[i+lda*j]);
		}
		printf(";\n");
	}
	printf("];\n");
}

void print_p_matrix(double* A, int lda, int n, int* iprm, char* name)
{
	int i, j;

	printf("%s=[", name);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%.18f, ", A[iprm[i]+lda*j]);
		}
		printf(";\n");
	}
	printf("];\n");
}

void print_p_lower_matrix(double* A, int lda, int n, int* iprm, char* name)
{
	int i, j;

	printf("%s=[", name);
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			printf("%.18f, ", A[iprm[i] + lda*j]);
		}
		printf("1, ");
		for (j = i + 1; j < n; j++) {
			printf("0, ");
		}
		printf(";\n");
	}
	printf("];\n");
}

void print_p_upper_matrix(double* A, int lda, int n, int* iprm, char* name)
{
	int i, j;

	printf("%s=[", name);
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			printf("0, ");
		}
		for (j = i; j < n; j++) {
			printf("%.18f, ", A[iprm[i] + lda*j]);
		}
		printf(";\n");
	}
	printf("];\n");
}

#define N	500
int n = 10;

double A[N*N];
double b[N], x[N];
int ipvt[N-1];
int iprm[N];

void test_LP()
{
	int k, info;

	for (k = 0; k < 20; k++) {
		rand_matrix(A, N, n);

		print_matrix(A, N, n, "A");

		info = dgefa(A, N, n, ipvt);
		if (info) {
			printf("Info=%d;\n", info);
			continue;
		}

		rand_vector(b, n);
		print_vector(b, n, "b");

		if (k%2) {
			dgesl(A, N, n, ipvt, b);
			print_vector(b, n, "x");
			printf("\nA*x-b\n\n");
		}
		else {
			dgeslt(A, N, n, ipvt, b);
			print_vector(b, n, "x");
			printf("\nA'*x-b\n\n");
		}
	}
}

void test_CYP()
{
	int k, info;

	for (k = 0; k < 20; k++) {
		rand_matrix(A, N, n);

		print_matrix(A, N, n, "A");

		info = CYP_dgefa(A, N, n, iprm);
		if (info) {
			printf("Info=%d;\n", info);
			continue;
		}

		rand_vector(b, n);
		print_vector(b, n, "b");

		if (k%2) {
			CYP_dgesl(A, N, n, iprm, b, x);
			print_vector(x, n, "x");
			printf("\nA*x-b\n\n");
		}
		else {
			CYP_dgeslt(A, N, n, iprm, b, x);
			print_vector(x, n, "x");
			printf("\nA'*x-b\n\n");
		}
	}
}

double time_rand_matrix(int n, int count)
{
	int k;
	time_t start_time, end_time;

	start_time = time(NULL);
	for (k = 0; k<count; k++) {
		rand_matrix(A, N, n);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

double time_LP_dgefa(int n, int count)
{
	int k;
	time_t start_time, end_time;
	
	start_time = time(NULL);
	for (k=0; k<count; k++) {
		rand_matrix(A, N, n);
		dgefa(A, N, n, ipvt);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

double time_CYP_dgefa(int n, int count)
{
	int k;
	time_t start_time, end_time;

	start_time = time(NULL);
	for (k = 0; k<count; k++) {
		rand_matrix(A, N, n);
		CYP_dgefa(A, N, n, iprm);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

double time_LP_dgesl(int n, int count)
{
	int k;
	time_t start_time, end_time;

	start_time = time(NULL);
	for (k = 0; k<count; k++) {
		rand_matrix(A, N, n);
		dgesl(A, N, n, ipvt, b);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

double time_CYP_dgesl(int n, int count)
{
	int k;
	time_t start_time, end_time;

	start_time = time(NULL);
	for (k = 0; k<count; k++) {
		rand_matrix(A, N, n);
		CYP_dgesl(A, N, n, iprm, b, x);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

double time_LP_dgeslt(int n, int count)
{
	int k;
	time_t start_time, end_time;

	start_time = time(NULL);
	for (k = 0; k<count; k++) {
		rand_matrix(A, N, n);
		dgeslt(A, N, n, ipvt, b);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

double time_CYP_dgeslt(int n, int count)
{
	int k;
	time_t start_time, end_time;

	start_time = time(NULL);
	for (k = 0; k<count; k++) {
		rand_matrix(A, N, n);
		CYP_dgeslt(A, N, n, iprm, b, x);
	}
	end_time = time(NULL);
	return (double)(end_time-start_time)/(double)count;
}

void compare_LP_CYP_times()
{
	int n, count;
	double f, time_rm;
	double time_lp_dgefa, time_cyp_dgefa;
	double time_lp_dgesl, time_cyp_dgesl, time_lp_dgeslt, time_cyp_dgeslt;

	for (n = 100; n<=N; n += 50) {
		f = (double)N/(double)n;
		count = (int)(f*f*400);
		time_rm = time_rand_matrix(n, count);
		time_lp_dgefa = time_LP_dgefa(n, count)-time_rm;
		time_cyp_dgefa = time_CYP_dgefa(n, count)-time_rm;
		time_lp_dgesl = time_LP_dgesl(n, count)-time_rm;
		time_cyp_dgesl = time_CYP_dgesl(n, count)-time_rm;
		time_lp_dgeslt = time_LP_dgeslt(n, count)-time_rm;
		time_cyp_dgeslt = time_CYP_dgeslt(n, count)-time_rm;
		printf("n=%d: time_rm=%f, time_lp_dgefa=%f, time_cyp_dgefa=%f, "
			"time_lp_dgesl=%f, time_cyp_dgesl=%f, time_lp_dgeslt=%f, time_cyp_dgeslt=%f\n",
			n, time_rm, time_lp_dgefa, time_cyp_dgefa, time_lp_dgesl, time_cyp_dgesl,
			time_lp_dgeslt, time_cyp_dgeslt);
	}
}

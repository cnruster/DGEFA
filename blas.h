// DGEFA+DGESL functions for Gaussian elimination linear equation solver
// Translated into C by Yiping Cheng, Beijing Jiaotong University, China
// Email: ypcheng@bjtu.edu.cn

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int idamax(int n, double const dx[], int incx);
void dscal(int n, double da, double dx[], int incx);
void daxpy(int n, double da, double const dx[], int incx, double dy[], int incy);
double ddot(int n, double const dx[], int incx, double const dy[], int incy);

int dgefa(double A[], int lda, int n, int ipvt[]);
void dgesl(double const A[], int lda, int n, int const ipvt[], double b[]);
void dgeslt(double const A[], int lda, int n, int const ipvt[], double b[]);


#ifdef __cplusplus
}
#endif

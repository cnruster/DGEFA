#pragma once

#ifdef __cplusplus
extern "C" {
#endif 

int idamax2(int n, double const col[], int const iprm[], int rstart);
int CYP_dgefa(double A[], int lda, int n, int iprm[]);
void CYP_dgesl(double const A[], int lda, int n, int const iprm[], double const b[], double x[]);
void CYP_dgeslt(double const A[], int lda, int n, int const iprm[], double b[], double x[]);

#ifdef __cplusplus
}
#endif 

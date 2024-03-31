#include <math.h>
#include <assert.h>
#include "blas_cyp.h"

// This version of IDAMAX is called by CYP_DGEFA
// written by Yiping Cheng, Beijing Jiaotong University
int idamax2(int n, double const col[], int const iprm[], int rstart)
{
    assert(n > 0);
    assert(0 <= rstart && rstart < n);

    double xmag, dmax = 0;
    int i, idm = rstart-1;

    for (i = rstart; i < n; i++) {
        if ((xmag = fabs(col[iprm[i]])) > dmax) {
            idm = i;
            dmax = xmag;
        }
    }

    return idm;
}


#define Aentry(i,j)        A[iprm[i]+lda*(j)]

int CYP_dgefa(double A[], int lda, int n, int iprm[])
{
    double pivot, * colk, * colj;
    int i, j, k, kpp, m, pi, pk, info;

    assert(n > 0 && lda >= n);

    for (k = 0; k < n; k++) {
        iprm[k] = k;
    }

    info = 0;

    // GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
    for (k = 0, kpp = 1; kpp < n; k = kpp++) {
        colk = &A[lda * k];
        // FIND m = PIVOT INDEX
        m = idamax2(n, colk, iprm, k);
        if (m < k) {
            // ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
            info = kpp;
            continue;
        }

        pk = iprm[m];
        if (m != k) { // INTERCHANGE IF NECESSARY
            iprm[m] = iprm[k];
            iprm[k] = pk;
        }

        // ROW ELIMINATION WITH COLUMN INDEXING
        pivot = colk[pk];
        for (i = kpp; i < n; i++) {
            colk[iprm[i]] /= pivot;
        }

        for (j = kpp; j < n; j++) {
            colj = &A[lda * j];
            for (i = kpp; i < n; i++) {
                pi = iprm[i];
                // Aentry(i,j) -= Aentry(i,k)*Aentry(k, j)
                colj[pi] -= colk[pi] * colj[pk];
            }
        }
    }

    if (Aentry(n-1, n-1) == 0.0) {
        info = n;
    }

    return info;
}

/*
CYP_DGESL solves the double precision system
A * X = B or TRANS(A)* X = B
using the factors computed by CYP_DGEFA.

IPRM    INTEGER(N)
the permutation formed during the Gaussian elimination process
*/

void
CYP_dgesl(double const A[], int lda, int n, int const iprm[], double const b[], double x[])
{
    double tmp;
    int k, pk, j;

    // SOLVE A*X=B, where A has been LU-decomposed

    // FIRST SOLVE L*Y = B. x[] stores Y
    for (k = 0; k<n; k++) {
        pk = iprm[k];
        tmp = b[pk];
        for (j = 0; j<k; j++) {
            tmp -= A[pk + lda*j] * x[j];
        }
        x[k] = tmp;
    }

    //  NOW SOLVE  U*X = Y
    for (k = n-1; k>=0; k--) {
        pk = iprm[k];
        tmp = x[k];
        for (j = k+1; j<n; j++) {
            tmp -= A[pk + lda*j] * x[j];
        }
        x[k] = tmp / A[pk + lda*k];
    }
}


/*
    CYP_dgeslt solves the double precision system
    A'*x = b
    using the factors computed by CYP_dgefa.
*/
void
CYP_dgeslt(double const A[], int lda, int n, int const iprm[], double b[], double x[])
{
    double tmp;
    double const* colk;
    int k, i, pi;

    // SOLVE trans(A)*X=B, where A has been LU-decomposed

    // FIRST SOLVE trans(U)*Y = B. b[] stores Y
    for (k = 0; k<n; k++) {
        tmp = b[k];
        colk = &A[lda*k];
        for (i = 0; i<k; i++) {
            tmp -= colk[iprm[i]]*b[i];
        }
        b[k] = tmp/A[iprm[k] + lda*k];
    }

    //  NOW SOLVE  trans(L)*X = Y
    for (k = n-1; k>=0; k--) {
        tmp = b[k];
        colk = &A[lda*k];
        for (i = k+1; i<n; i++) {
            pi = iprm[i];
            tmp -= colk[pi]*x[pi];
        }
        x[iprm[k]] = tmp;
    }
}
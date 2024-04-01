// DGEFA+DGESL functions for Gaussian elimination linear equation solver
// Translated into C by Yiping Cheng, Beijing Jiaotong University, China
// Email: ypcheng@bjtu.edu.cn

#include <math.h>
#include <assert.h>
#include "blas.h"

/***BEGIN PROLOGUE  IDAMAX
***PURPOSE  Find the smallest index of that component of a vector
having the maximum magnitude.
***CATEGORY  D1A2
***TYPE      DOUBLE PRECISION(ISAMAX - S, IDAMAX - D, ICAMAX - C)
***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
***AUTHOR  Lawson, C.L., (JPL)
Hanson, R.J., (SNLA)
Kincaid, D.R., (U.of Texas)
Krogh, F.T., (JPL)
***DESCRIPTION

B L A S  Subprogram
Description of Parameters

--Input--
N  number of elements in input vector(s)
DX  double precision vector with N elements
INCX  storage spacing between elements of DX

--Output--
IDAMAX  smallest index(zero if N.LE. 0)

Find smallest index of maximum magnitude of double precision DX.
IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX + (I - 1) * INCX)),
where IX = 1 if INCX.GE. 0, else IX = 1 + (1 - N) * INCX.

***REFERENCES  C.L.Lawson, R.J.Hanson, D.R.Kincaidand F.T.
Krogh, Basic linear algebra subprograms for Fortran
usage, Algorithm No. 539, Transactions on Mathematical
Software 5, 3 (September 1979), pp. 308 - 323.
***ROUTINES CALLED(NONE)
***REVISION HISTORY(YYMMDD)
791001  DATE WRITTEN
890531  Changed all specific intrinsics to generic. (WRB)
890531  REVISION DATE from Version 3.2
891214  Prologue converted to Version 4.0 format. (BAB)
900821  Modified to correct problem with a negative increment.
920501  Reformatted the REFERENCES section. (WRB)
240327  Translated into C by Yiping Cheng. (CYP)
***END PROLOGUE  IDAMAX
*/
int idamax(int n, double const dx[], int incx)
{
    if (n <= 1) {
        return n-1;
    }

    double dmax, xmag;
    int i, ix;
    int idm = 0;

    if (incx != 1) {
        // Code for increments not equal to 1.
        ix = (incx < 0) ? (1 - n) * incx : 0;
        dmax = fabs(dx[ix]);

        for (i = 1; i < n; i++) {
            ix += incx;
            if ((xmag = fabs(dx[ix])) > dmax) {
                idm = i;
                dmax = xmag;
            }
        }
    }
    else {
        // Code for increments equal to 1.
        dmax = fabs(dx[0]);
        for (i = 1; i < n; i++) {
            if ((xmag = fabs(dx[i])) > dmax) {
                idm = i;
                dmax = xmag;
            }
        }
    }

    return idm;
}

/***BEGIN PROLOGUE  DSCAL
***PURPOSE  Multiply a vector by a constant.
***CATEGORY  D1A6
***TYPE      DOUBLE PRECISION(SSCAL - S, DSCAL - D, CSCAL - C)
***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
***AUTHOR  Lawson, C.L., (JPL)
Hanson, R.J., (SNLA)
Kincaid, D.R., (U.of Texas)
Krogh, F.T., (JPL)
***DESCRIPTION

B L A S  Subprogram
Description of Parameters

--Input--
N  number of elements in input vector(s)
DA  double precision scale factor
DX  double precision vector with N elements
INCX  storage spacing between elements of DX

--Output--
DX  double precision result(unchanged if N.LE.0)

Replace double precision DX by double precision DA*DX.
For I = 0 to N - 1, replace DX(IX + I * INCX) with  DA* DX(IX + I * INCX),
where IX = 1 if INCX.GE. 0, else IX = 1 + (1 - N) * INCX.

***REFERENCES  C.L.Lawson, R.J.Hanson, D.R.Kincaidand F.T.
Krogh, Basic linear algebra subprograms for Fortran
usage, Algorithm No. 539, Transactions on Mathematical
Software 5, 3 (September 1979), pp. 308 - 323.
***ROUTINES CALLED(NONE)
***REVISION HISTORY(YYMMDD)
791001  DATE WRITTEN
890831  Modified array declarations. (WRB)
890831  REVISION DATE from Version 3.2
891214  Prologue converted to Version 4.0 format. (BAB)
900821  Modified to correct problem with a negative increment.(WRB)
920501  Reformatted the REFERENCES section. (WRB)
240327  Translated into C by Yiping Cheng. (CYP)
***END PROLOGUE  DSCAL
*/
void dscal(int n, double da, double dx[], int incx)
{
    if (n <= 0) {
    }
    else if (incx != 1) {
        // Code for increment not equal to 1.
        int i, ix;
        ix = (incx < 0) ? (1 - n) * incx : 0;
        for (i = 0; i < n; i++) {
            dx[ix] *= da;
            ix += incx;
        }
    }
    else {
        int i, m = n % 5;
        for (i = 0; i < m; i++) {
            dx[i] *= da;
        }
        // Clean-up loop so remaining vector length is a multiple of 5.
        for (; i < n; i += 5) {
            dx[i] *= da;
            dx[i+1] *= da;
            dx[i+2] *= da;
            dx[i+3] *= da;
            dx[i+4] *= da;
        }
    }
}

/*
***BEGIN PROLOGUE  DDOT
***PURPOSE  Compute the inner product of two vectors.
***CATEGORY  D1A4
***TYPE      DOUBLE PRECISION(SDOT - S, DDOT - D, CDOTU - C)
***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
***AUTHOR  Lawson, C.L., (JPL)
Hanson, R.J., (SNLA)
Kincaid, D.R., (U.of Texas)
Krogh, F.T., (JPL)
***DESCRIPTION

B L A S  Subprogram
Description of Parameters

--Input--
N  number of elements in input vector(s)
DX  double precision vector with N elements
INCX  storage spacing between elements of DX
DY  double precision vector with N elements
INCY  storage spacing between elements of DY

--Output--
DDOT  double precision dot product(zero if N.LE. 0)

Returns the dot product of double precision DX and DY.
DDOT = sum for I = 0 to N - 1 of  DX(LX + I * INCX) * DY(LY + I * INCY),
where LX = 1 if INCX.GE. 0, else LX = 1 + (1 - N) * INCX, and LY is
defined in a similar way using INCY.

***REFERENCES  C.L.Lawson, R.J.Hanson, D.R.Kincaidand F.T.
Krogh, Basic linear algebra subprograms for Fortran
usage, Algorithm No. 539, Transactions on Mathematical
Software 5, 3 (September 1979), pp. 308 - 323.
***ROUTINES CALLED(NONE)
***REVISION HISTORY(YYMMDD)
791001  DATE WRITTEN
890831  Modified array declarations.  (WRB)
890831  REVISION DATE from Version 3.2
891214  Prologue converted to Version 4.0 format.  (BAB)
920310  Corrected definition of LX in DESCRIPTION.  (WRB)
920501  Reformatted the REFERENCES section.  (WRB)
240327  Translated into C by Yiping Cheng.  (CYP)
***END PROLOGUE  DDOT
*/
double ddot(int n, double const dx[], int incx, double const dy[], int incy)
{
    double dotprod = 0.0;

    if (n <= 0) {
    }
    else if (incx != incy || incx <= 0) {
        // Code for unequal or nonpositive increments
        int ix, iy, i;

        ix = (incx < 0) ? (1 - n) * incx : 0;
        iy = (incy < 0) ? (1 - n) * incy : 0;

        for (i = 0; i < n; i++) {
            dotprod += dx[ix] * dy[iy];
            ix += incx;
            iy += incy;
        }
    }
    else if (incx == 1) {
        // Code for both increments equal to 1
        int i, m = n % 5;
        for (i = 0; i < m; i++) {
            dotprod += dx[i] * dy[i];
        }

        // Clean-up loop so remaining vector length is a multiple of 5
        for (; i < n; i += 5) {
            dotprod += dx[i] * dy[i] + dx[i+1] * dy[i+1] + dx[i+2] * dy[i+2]
                + dx[i+3] * dy[i+3] + dx[i+4] * dy[i+4];
        }
    }
    else {
        // Code for equal, positive, non-unit increments
        int i, ns = n * incx;
        for (i = 0; i < ns; i += incx) {
            dotprod += dx[i] * dy[i];
        }
    }

    return dotprod;
}

/***BEGIN PROLOGUE  DAXPY
***PURPOSE  Compute a constant times a vector plus a vector.
***CATEGORY  D1A7
***TYPE      DOUBLE PRECISION(SAXPY - S, DAXPY - D, CAXPY - C)
***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
***AUTHOR  Lawson, C.L., (JPL)
Hanson, R.J., (SNLA)
Kincaid, D.R., (U.of Texas)
Krogh, F.T., (JPL)
***DESCRIPTION

B L A S  Subprogram
Description of Parameters

--Input--
N  number of elements in input vector(s)
DA  double precision scalar multiplier
DX  double precision vector with N elements
INCX  storage spacing between elements of DX
DY  double precision vector with N elements
INCY  storage spacing between elements of DY

--Output--
DY  double precision result(unchanged if N.LE. 0)

Overwrite double precision DY with double precision DA* DX + DY.
For I = 0 to N - 1, replace  DY(LY + I * INCY) with DA * DX(LX + I * INCX) +
DY(LY + I * INCY),
where LX = 1 if INCX.GE. 0, else LX = 1 + (1 - N) * INCX, and LY is
defined in a similar way using INCY.

***REFERENCES  C.L.Lawson, R.J.Hanson, D.R.Kincaidand F.T.
Krogh, Basic linear algebra subprograms for Fortran
usage, Algorithm No. 539, Transactions on Mathematical
Software 5, 3 (September 1979), pp. 308 - 323.
***ROUTINES CALLED(NONE)
***REVISION HISTORY(YYMMDD)
791001  DATE WRITTEN
890831  Modified array declarations.  (WRB)
890831  REVISION DATE from Version 3.2
891214  Prologue converted to Version 4.0 format.  (BAB)
920310  Corrected definition of LX in DESCRIPTION.  (WRB)
920501  Reformatted the REFERENCES section.  (WRB)
240327  Translated into C by Yiping Cheng.  (CYP)
***END PROLOGUE  DAXPY
*/
void daxpy(int n, double da, double const dx[], int incx, double dy[], int incy)
{
    if (n <= 0 || da == 0.0) {
    }
    else if (incx != incy || incx <= 0) {
        // Code for unequal or nonpositive increments
        int ix, iy, i;

        ix = (incx < 0) ? (1 - n) * incx : 0;
        iy = (incy < 0) ? (1 - n) * incy : 0;

        for (i = 0; i < n; i++) {
            dy[iy] += da*dx[ix];
            ix += incx;
            iy += incy;
        }
    }
    else if (incx == 1) {
        // Code for both increments equal to 1
        int i, m = n % 4;
        for (i = 0; i < m; i++) {
            dy[i] += da*dx[i];
        }

        // Clean-up loop so remaining vector length is a multiple of 4
        for (; i < n; i += 4) {
            dy[i] += da*dx[i];
            dy[i+1] += da*dx[i+1];
            dy[i+2] += da*dx[i+2];
            dy[i+3] += da*dx[i+3];
        }
    }
    else {
        // Code for equal, positive, non-unit increments.
        int i, ns = n * incx;
        for (i = 0; i < ns; i += incx) {
            dy[i] += da*dx[i];
        }
    }
}

/*
dgefa factors a double precision matrix by gaussian elimination.

    on entry

        A       double precision(lda, n)
                the matrix to be factored.

        lda     integer
                the leading dimension of the array A

        n       integer
                the order of the matrix A

    on return

        A       an upper triangular matrix and the multipliers
                which were used to obtain it.
                the factorization can be written  A = L*U  where
                L is a product of permutation and unit lower
                triangular matrices and U is upper triangular.

        ipvt    integer(n-1)
                an integer vector of pivot indices.

        info    integer
                = 0  normal value.
                = k  if  u(k,k) .eq. 0.0 . This is not an error
                    condition for this subroutine, but it does
                    indicate that dgesl or dgedi will divide by zero
                    if called.

    linpack. this version dated 08/14/78 .
    cleve moler, university of new mexico, argonne national lab.
    Translated into C by Yiping Cheng, 27/03/2024
*/

int dgefa(double A[], int lda, int n, int ipvt[])
{
    double tmp, *colk, *colj;
    int k, kpp, m, j;
    int info = 0;

    // Gaussian elimination with partial pivoting

    for (k=0, kpp=1; kpp < n; k=kpp++) {
        colk = &A[lda*k];
        // Find pivot index
        m = k + idamax(n-k, &colk[k], 1);
        ipvt[k] = m;

        if ((tmp = colk[m]) == 0.) {
            // Zero pivot implies this row already triangularized
            info = kpp;
            continue;
        }

        if (m != k) {
            // Interchange
            colk[m] = colk[k];
            colk[k] = tmp;

            // Compute multipliers and scale
            dscal(n-kpp, -1./tmp, &colk[kpp], 1);

            // Column elimination with row indexing
            for (j = kpp; j<n; j++) {
                colj = &A[lda*j];
                tmp = colj[m];
                colj[m] = colj[k];
                colj[k] = tmp;
                daxpy(n-kpp, tmp, &colk[kpp], 1, &colj[kpp], 1);
            }
        }
        else {
            // Compute multipliers and scale
            dscal(n-kpp, -1./tmp, &colk[kpp], 1);

            // Column elimination with row indexing
            for (j = kpp; j<n; j++) {
                colj = &A[lda*j];
                daxpy(n-kpp, colj[m], &colk[kpp], 1, &colj[kpp], 1);
            }
        }
    }

    // now k becomes n-1
    if (A[k + lda*k] == 0.)
        info = n;

    return info;
}

/*
    dgesl solves the double precision system
    A*x = b
    using the factors computed by dgefa.

    on entry

        A       double precision(lda, n)
                the output from dgefa.

        n       integer
                the order of the matrix  A .

        ipvt    integer(n-1)
                the pivot vector from dgefa.

        b       double precision(n)
                the right hand side vector.

    on return

        b       the solution vector x.

    error condition

        A division by zero will occur if the input factor contains A
        zero on the diagonal. technically this indicates singularity
        but it is often caused by improper arguments or improper
        setting of lda. it will not occur if the subroutines are
        called correctly and if dgefa has set info.eq. 0 .

    linpack. this version dated 08/14/78 .
    cleve moler, university of new mexico, argonne national lab.
    Translated into C by Yiping Cheng, 27/03/2024
*/

void dgesl(double const A[], int lda, int n, int const ipvt[], double b[])
{
    double tmp;
    double const *colk;
    int k, kpp, m;

    // solve  A * x = b

    // first solve  L*y = b
    for (k = 0, kpp = 1; kpp < n; k = kpp++) {
        m = ipvt[k];
        tmp = b[m];
        if (m!=k) {
            b[m] = b[k];
            b[k] = tmp;
        }

        daxpy(n-kpp, tmp, &A[kpp + lda*k], 1, &b[kpp], 1);
    }

    // now solve  U*x = y
    for (; k>=0; k--) {
        colk = &A[lda*k];
        b[k] /= colk[k];
        daxpy(k, -b[k], colk, 1, b, 1);
    }
}


/*
    dgeslt solves the double precision system
    A'*x = b
    using the factors computed by dgefa.
    Translated into C by Yiping Cheng, 27/03/2024
*/
void dgeslt(double const A[], int lda, int n, int const ipvt[], double b[])
{
    int k, kpp, m;
    double tmp;
    double const *colk;

    // solve trans(A)*x = b

    //  First solve trans(U)*y = b
    for (k=0; k<n; k++) {
        colk = &A[lda*k];
        tmp = ddot(k, colk, 1, b, 1);
        b[k] = (b[k] - tmp) / colk[k];
    }

    // Now solve trans(L)* x = y
    for (kpp=n-1, k=n-2; k>=0; kpp=k--) {
        b[k] += ddot(n-kpp, &A[kpp + lda*k], 1, &b[kpp], 1);
        m = ipvt[k];
        if (m != k) {
            tmp = b[m];
            b[m] = b[k];
            b[k] = tmp;
        }
    }
}
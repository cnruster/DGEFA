**DGEFA : Gaussian elimination (LU decomposition) with partial pivoting**

**DGESL and DGESLT : Backsubstitution**

Together they can solve any general linear system.


**This library is the highest quality code when you need to solve a general linear system!**

**You get a deep understanding of Gaussian elimination by reading this code**

It is even (slightly) better than the orginal LinPack code which is in Fortran.

**Cautions for using DGEFA+DGESL+DGESLT into your application**
1) The functions are translated from Fortran, so the matrices should be stored in column-major mode. That is, the address of A[i,j] is &A[i+lda*j]
2) The indices start with 0, not 1 as in Fortran.
3) lda is the leading dimension of a matrix A. Usually lda is the number of rows of A, but not always. If A is stored as a submatrix of a matrix B, then lda may be the number of rows of B.
4) Dgefa destroys the original matrix A. Upon return, the upper triangular+diagonal part of A becomes U, and the lower triangular part of A contains the Gaussian elimination coefficients (but not L literally).
5) info is now a return value of dgefa, not an argument. info==0 means A is nonsingular, and info!=0 means A is singular.
6) When info==0, do not call dgesl or dgeslt! The linear system is ill-conditioned. And a divide-by-zero exception will occur if you insist on calling them.
7) I have split the original DGESL subroutine into two functions: dgesl and dgeslt. So we now don't need the job argument.
8) My numerical experiment shows CYP_dgefa is less efficient than dgefa. I have not yet found the reason for this, as I originally thought CYP_dgefa would outperform dgefa with much fewer swappings. Anyway, it is now clear that the standard dgefa must be preferred.

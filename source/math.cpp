#pragma once


#include <cstdio>
#include <cstdlib>
#include <cstring>


// BLAS / LAPACK declarations

extern "C"
{
  void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
  void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
  void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
}


// Inverts n x n matrix a using LAPACK. Returns the result in inv.

void invert_matrix(double* a, int n, double* inv)
{
  int*    ipiv = new int[n];
  double* work = new double[n];  // minimum required size

  memcpy(inv, a, n * n * sizeof(double));

  int info;
  dgetrf_(&n, &n, inv, &n, ipiv, &info);
  if (info != 0) { printf("Matrix inversion failed!\n"); exit(1); }
  dgetri_(&n, inv, &n, ipiv, work, &n, &info);
  if (info != 0) { printf("Matrix inversion failed!\n"); exit(1); }

  delete[] ipiv;
  delete[] work;
}


// Multiplies matrix a of size m x k with matrix b of size k x n. The result is
// returned in m x n matrix c.
//
// To use the column-major routine directly, this routine uses the fact that
// AB = C transposes to B^T A^T = C^T (the output is column major and thus
// always transposed, so the input needs to change to compensate), and A and B
// are already viewed as transposed due to row major ordering.

void matmul(int m, int n, int k, double alpha, double* a, double beta,
            double* b, double* c)
{
  char transp = 'n';

  dgemm_(&transp, &transp, &n, &m, &k, &alpha, b, &n, a, &k, &beta, c, &n);
}

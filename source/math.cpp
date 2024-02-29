#pragma once


#include <cstdio>
#include <cstdlib>
#include <cstring>

// BLAS / LAPACK declarations

extern "C"
{
  void dgetrf_(int m, int n, double* a, int lda, int* ipiv, int info);
  void dgetri_(int n, double* a, int lda, int* ipiv, double* work, int lwork, int info);
}


// Inverts n x n matrix a using LAPACK. Returns the result in inv.
void invert_matrix(double* a, int n, double* inv)
{
  int*    ipiv = new int[n];
  double* work = new double[n];  // minimum required size

  memcpy(inv, a, n * n * sizeof(double));

  int info;
  dgetrf_(n, n, inv, n, ipiv, info);
  if (info != 0) { printf("Matrix inversion failed!\n"); exit(1); }
  dgetri_(n, inv, n, ipiv, work, n, info);
  if (info != 0) { printf("Matrix inversion failed!\n"); exit(1); }

  delete[] ipiv;
  delete[] work;
}

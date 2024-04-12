#pragma once


#include "linalg.cpp"
#include "linear_solver.cpp"


struct test_system
{
  double* A;
  usize   m;

  void matvec(double* v, double* mv);

  // --

  test_system(usize w, usize h);
  ~test_system();

  test_system(const test_system& oth)            = delete;
  test_system(test_system&& oth)                 = delete;
  test_system& operator=(const test_system& oth) = delete;
  test_system& operator=(test_system&& oth)      = delete;
};


void test_system::matvec(double* v, double* mv)
{
  matmul(m, 1, m, 1., A, v, 0., mv);
}


test_system::test_system(usize w, usize h)
{
  m = w * h;
  A = new double[m * m];

  memset(A, 0, m * m * sizeof(double));

  usize Arow = 0;
  for (usize hi = 0; hi < h; ++hi)
  {
    for (usize wi = 0; wi < w; ++wi)
    {
      A[m * (Arow) + (Arow)] = +4.;

      if (wi > 0)
        A[m * (Arow) + (Arow - 1)] = -1.;
      if (wi + 1 < w)
        A[m * (Arow) + (Arow + 1)] = -1.;
      if (hi > 0)
        A[m * (Arow) + (Arow - w)] = -1.;
      if (hi + 1 < h)
        A[m * (Arow) + (Arow + w)] = -1.;

      ++Arow;
    }
  }
}


test_system::~test_system()
{
  delete[] A;
}


int test_GMRES()
{
  usize w = 10;
  usize h = 10;

  test_system sys(10, 10);

  double* x   = new double[sys.m];
  double* b   = new double[sys.m];
  double* tmp = new double[sys.m];

  // manufacture solution to 1
  for (usize i = 0; i < sys.m; ++i) { x[i] = 1.; }
  sys.matvec(x, b);

  // testing...

  memset(x, 0, sys.m * sizeof(double));

  GMRES_settings set;
  set.ninner = 50;
  set.nouter = 10;
  set.rtol   = 1e-14;

  GMRES_workspace wsp(set, sys.m);

  GMRES(sys, b, x, sys.m, set, wsp);

  // test linear residual
  sys.matvec(x, tmp);
  for (usize i = 0; i < sys.m; ++i) { tmp[i] -= b[i]; }
  double real_resid = l2norm(tmp, sys.m);
  printf("real resid: %+.17e\n", real_resid);


  bool test_fail = false;
  for (usize i = 0; i < sys.m; ++i)
  {
    double err = std::abs(1. - x[i]);
    if (err > 1e-14)
    {
      test_fail = true;
      printf("Solution deviation from 1. of %+17e!\n", err);
    }
  }

  delete[] x;
  delete[] b;
  delete[] tmp;

  return (test_fail) ? 1 : 0;
}

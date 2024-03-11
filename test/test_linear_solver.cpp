#pragma once


#include "math.cpp"
#include "linear_solver.cpp"


usize   m__;
double* A__;


void matvec(double* v, double* res)
{
  matmul(m__, 1, m__, 1., A__, v, 0., res);
}


int test_GMRES()
{
  usize w = 10;
  usize h = 10;

  m__ = w * h;

          A__ = new double[m__ * m__];
  double* x   = new double[m__];
  double* b   = new double[m__];

  // system setup

  usize Arow = 0;
  for (usize hi = 0; hi < h; ++hi)
  {
    for (usize wi = 0; wi < w; ++wi)
    {
      A__[m__ * (Arow) + (Arow)] = +4.;

      if (wi > 0)
        A__[m__ * (Arow) + (Arow - 1)] = -1.;
      if (wi + 1 < w)
        A__[m__ * (Arow) + (Arow + 1)] = -1.;
      if (hi > 0)
        A__[m__ * (Arow) + (Arow - w)] = -1.;
      if (hi + 1 < h)
        A__[m__ * (Arow) + (Arow + w)] = -1.;

      ++Arow;
    }
  }

  for (usize i = 0; i < m__; ++i) { x[i] = 1.; }
  matmul(m__, 1, m__, 1., A__, x, 0., b);

  // testing...

  memset(x, 0, m__ * sizeof(double));

  GMRES(matvec, b, m__, x, 100, 1, 1e-10);

  delete[] A__;
  delete[] x;
  delete[] b;

  return 0;
}

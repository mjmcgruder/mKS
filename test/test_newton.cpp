#pragma once


#include <cfloat>

#include "time.cpp"


struct bowl
{
  // external facing vars

  double F[2];

  // worksapce and references

  double* xk;     // linearization state (not owned)
  double  Fp[2];  // perturbed residual
  double  xp[2];  // perturbed state

  // --

  void uresid(double* u, double* F);
  void matvec(double* v, double* mv);
  void update_guess(double* x);
};


void bowl::update_guess(double* u)
{
  xk = u;
  uresid(u, F);
}


void bowl::uresid(double* u, double* F)
{
  double x   = u[0];
  double y   = u[1];
  double xm1 = x - 1.;
  double ym1 = y - 1.;
  F[0] = xm1 * xm1 + ym1 * ym1 - 1.;
  F[1] = x - y;
}


void bowl::matvec(double* v, double* mv)
{
  double eps  = 1e-8;
  double ieps = 1. / eps;

  for (usize i = 0; i < 2; ++i)
  {
    xp[i] = xk[i] + eps * v[i];
  }

  uresid(xp, Fp);

  for (usize i = 0; i < 2; ++i)
  {
    mv[i] = (Fp[i] - F[i]) * ieps;
  }
}


int test_newton()
{
  double tol = 1e-10;
  bowl system;

  double x[2] = {2., 1.};

  newton_settings  nset{2, 100, 1e-10, 1e-10, 10};
  newton_workspace nwrk(2);

  newton_solve(system, x, 2, nset, nwrk);

  printf("Answer: %+.17e %+.17e\n", x[0], x[1]);

  double err0 = x[0] - (1. + 0.5 * sqrt(2));
  double err1 = x[1] - (1. + 0.5 * sqrt(2));

  printf("error: %+.17e %+.17e\n", err0, err1);
  if (err0 > tol || err1 > tol) { return 1; }

  return 0;
}

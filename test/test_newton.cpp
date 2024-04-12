#pragma once


#include <cfloat>

#include "time.cpp"
#include "time_schemes.cpp"


struct bowl
{
  void uresid(double* u, double* F);
};


void bowl::uresid(double* u, double* F)
{
  double x   = u[0];
  double y   = u[1];
  double xm1 = x - 1.;
  double ym1 = y - 1.;
  F[0] = xm1 * xm1 + ym1 * ym1 - 1.;
  F[1] = x - y;
}


int test_newton()
{
  double tol = 1e-10;
  double x[2] = {2., 1.};

  newton_settings  nset{};
  nset.gmres.ninner = 2;
  nset.gmres.nouter = 100;
  nset.gmres.rtol   = 1e-8;
  nset.gmres.ctol   = 1e-11;
  nset.nonlin_tol   = 1e-10;
  nset.max_iter     = 10;

  newton_workspace nwrk(2);
  GMRES_workspace gwrk(nset.gmres, 2);

  bowl system;
  matvec_fd<bowl> linsys(system, 2);

  newton_solve(linsys, x, 2, nset, gwrk, nwrk);

  printf("Answer: %+.17e %+.17e\n", x[0], x[1]);

  double err0 = x[0] - (1. + 0.5 * sqrt(2));
  double err1 = x[1] - (1. + 0.5 * sqrt(2));

  printf("error: %+.17e %+.17e\n", err0, err1);
  if (err0 > tol || err1 > tol) { return 1; }

  return 0;
}

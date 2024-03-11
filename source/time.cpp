#pragma once


#include <cmath>

#include "linalg.cpp"
#include "residual.cpp"
#include "linear_solver.cpp"


void euler(comp_geom& geom, usize nt, double dt, double* state)
{
  double* resid = new double[geom.size_st()];
  double* f     = new double[geom.size_st()];

  for (usize it = 0; it < nt; ++it)
  {
    compute_residual(geom, state, resid, f);

    printf("%zu: %+.3e\n", it, norm_perelem(f, geom.size_st()));

    for (usize i = 0; i < geom.size_st(); ++i)
    {
      state[i] += dt * f[i];
    }
  }

  delete[] resid;
  delete[] f;
}


void rk4(comp_geom& geom, usize nt, double dt, double* y1, double* hist)
{
  double* resid = new double[geom.size_st()];

  double* k1 = new double[geom.size_st()];
  double* k2 = new double[geom.size_st()];
  double* k3 = new double[geom.size_st()];
  double* k4 = new double[geom.size_st()];

  double* y2 = new double[geom.size_st()];
  double* y3 = new double[geom.size_st()];
  double* y4 = new double[geom.size_st()];

  for (usize it = 0; it < nt; ++it)
  {
    memcpy(hist + it * geom.size_st(), y1, geom.size_st() * sizeof(double));

    compute_residual(geom, y1, resid, k1);

    printf("%zu: %+.3e\n", it, norm_perelem(k1, geom.size_st()));

    for (usize i = 0; i < geom.size_st(); ++i) y2[i] = y1[i] + 0.5 * dt * k1[i];

    compute_residual(geom, y2, resid, k2);

    for (usize i = 0; i < geom.size_st(); ++i) y3[i] = y1[i] + 0.5 * dt * k2[i];

    compute_residual(geom, y3, resid, k3);

    for (usize i = 0; i < geom.size_st(); ++i) y4[i] = y1[i] + dt * k3[i];

    compute_residual(geom, y4, resid, k4);

    for (usize i = 0; i < geom.size_st(); ++i)
    {
      y1[i] += (dt / 6.) * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    }
  }

  memcpy(hist + nt * geom.size_st(), y1, geom.size_st() * sizeof(double));

  delete[] resid;

  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] k4;

  delete[] y2;
  delete[] y3;
  delete[] y4;
}


void BDF1(double* unp1, double* un, double dt, comp_geom& geom, double* resid,
          double* f, double* F)
{
  // compute steady residual
  compute_residual(geom, unp1, resid, f);

  // apply time scheme
  for (usize i = 0; i < geom.size_st(); ++i)
  {
    F[i] = unp1[i] - un[i] - dt * f[i];
  }
}


double* NEWTON_uk;  // current nonlinear solve guess


// BDF1 callback supporting globals
double*    BDF1_un;       // known (previous) state
double*    BDF1_uk_pert;  // perturbed nonlinear solve guess (filled internally)
double     BDF1_dt;       // }
comp_geom* BDF1_geom;     // | residual data / workspace
double*    BDF1_resid;    // |
double*    BDF1_f;        // }
double*    BDF1_Fo;       // original unsteady residual (assigned before entry)
double*    BDF1_Fp;       // perturbed unsteady residual


void BDF1_matvec(double* v, double* mv)
{
  double normv = norm(v, BDF1_geom->size_st());
  double eps   = 1e-8;

  // divide by norm for pure direction, multiply by eps for length

  for (usize i = 0; i < BDF1_geom->size_st(); ++i)
  {
    BDF1_uk_pert[i] = NEWTON_uk[i] + eps * (v[i] / normv);
  }

  BDF1(BDF1_uk_pert, BDF1_un, BDF1_dt, *BDF1_geom, BDF1_resid, BDF1_f, BDF1_Fp);

  for (usize i = 0; i < BDF1_geom->size_st(); ++i)
  {
    mv[i] = (BDF1_Fp[i] - BDF1_Fo[i]) / eps;
  }
}


void newton(comp_geom& geom, double dt, double* un)
{
  BDF1_dt   = dt;
  BDF1_geom = &geom;

  NEWTON_uk  = new double[geom.size_st()];
  BDF1_resid = new double[geom.size_st()];
  BDF1_f     = new double[geom.size_st()];
  BDF1_Fo    = new double[geom.size_st()];
  BDF1_Fp    = new double[geom.size_st()];
  double* x  = new double[geom.size_st()];
  double* b  = new double[geom.size_st()];

  // BEGIN NONLINEAR SOLUTION ITERATION

  memcpy(NEWTON_uk, un, geom.size_st() * sizeof(double));
  memcpy(BDF1_un, un, geom.size_st() * sizeof(double));

  BDF1(NEWTON_uk, un, dt, geom, BDF1_resid, BDF1_f, BDF1_Fo);  // call function for initial RHS
  double nonlinear_residual = norm(BDF1_Fo, geom.size_st());

  while (nonlinear_residual > 1e-8)
  {
    for (usize i = 0; i < geom.size_st(); ++i) { b[i] = -BDF1_Fo[i]; }
    memset(x, 0, geom.size_st() * sizeof(double));

    // solve with GMRES
    GMRES(BDF1_matvec, b, geom.size_st(), x, 100, 10, 1e-8);

    // update guess with (full?) step
    for (usize i = 0; i < geom.size_st(); ++i) { NEWTON_uk[i] += x[i]; }

    // evaluate next RHS
    BDF1(NEWTON_uk, un, dt, geom, BDF1_resid, BDF1_f, BDF1_Fo);
    nonlinear_residual = norm(BDF1_Fo, geom.size_st());
  }

  // return guess as next state
  memcpy(un, NEWTON_uk, geom.size_st() * sizeof(double));

  // END NONLINEAR SOLUTION ITERATION

  delete[] NEWTON_uk;
  delete[] BDF1_resid;
  delete[] BDF1_f;
  delete[] BDF1_Fo;
  delete[] BDF1_Fp;
  delete[] x;
  delete[] b;
}

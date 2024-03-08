#pragma once


#include <cmath>

#include "residual.cpp"


double norm2(double* a, usize len)
{
  double sum = 0.;
  for (usize i = 0; i < len; ++i)
  {
    sum += a[i] * a[i];
  }

  return sqrt((1. / double(len)) * sum);
}


void euler(comp_geom& geom, usize nt, double dt, double* state)
{
  double* resid = new double[geom.size_st()];
  double* f     = new double[geom.size_st()];

  for (usize it = 0; it < nt; ++it)
  {
    compute_residual(geom, state, resid, f);

    printf("%zu: %+.3e\n", it, norm2(f, geom.size_st()));

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

    printf("%zu: %+.3e\n", it, norm2(k1, geom.size_st()));

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

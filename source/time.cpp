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


void rk4()
{
}

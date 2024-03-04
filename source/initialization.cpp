#pragma once


#include "geometry.cpp"


void init_constant(comp_geom& geom, double val, double* state)
{
  for (usize i = 0; i < geom.size_st(); ++i)
  {
    state[i] = val;
  }
}


void init_gaussian(comp_geom& geom, double a, double b, double c, double* state)
{
  double* n = new double[geom.p + 1];
  lagrange_nodes(geom.p, n);

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    for (usize ni = 0; ni < geom.p + 1; ++ni)
    {
      double x                   = geom.ref2glo(ei, n[ni]);
      double g                   = a * exp(-((x - b) * (x - b)) / (2. * c * c));
      geom.idx_st(state, ei, ni) = g;
    }
  }

  delete[] n;
}

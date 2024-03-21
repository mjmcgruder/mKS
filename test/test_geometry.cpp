#pragma once


#include <cstdio>

#include "types.cpp"
#include "geometry.cpp"
#include "linalg.cpp"


int test_mass_matrix()
{
  comp_geom geom(1., 3, 7, 1., 1., 1., 1., 1., 1.);

  double* M = new double[(geom.p + 1) * (geom.p + 1)];

  invert_matrix(geom.Minv, (geom.p + 1), M);

  double sum = 0.;
  for (usize i = 0; i < (geom.p + 1) * (geom.p + 1); ++i)
  {
    sum += M[i];
  }

  double err = std::abs(sum - 2.);
  if (err > 1e-14)
  {
    printf("Sum of mass matrix elements does not match length of reference "
           "domain! Got %+.17e, off by %+.17e.",
           sum, err);
    return 1.;
  }

  delete[] M;

  return 0;
}

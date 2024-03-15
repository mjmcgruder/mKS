#pragma once


#include <cmath>

#include "types.cpp"


double norm(double* a, usize len)
{
  double sum = 0.;
  for (usize i = 0; i < len; ++i)
  {
    sum += a[i] * a[i];
  }

  return sqrt(sum);
}

double norm_perelem(double* a, usize len)
{
  double sum = 0.;
  for (usize i = 0; i < len; ++i)
  {
    sum += a[i] * a[i];
  }

  return sqrt((1. / double(len)) * sum);
}


double dot(double* a, double* b, usize len)
{
  double sum = 0.;
  for (usize i = 0; i < len; ++i) { sum += a[i] * b[i]; }

  return sum;
}

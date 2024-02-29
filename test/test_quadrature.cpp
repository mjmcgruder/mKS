#pragma once


#include <cstdio>
#include <cmath>

#include "quadrature.cpp"


double polynomial(usize order, double x)
{
  double coeff[] = {0.49330289, 0.90463143, 0.18816676, 0.30646078, 0.18025477, 0.77162931, 0.03942296, 0.33686216, 0.26149018, 0.49742156, 0.47202936, 0.32424461, 0.25503326, 0.26378048, 0.50724305, 0.94544638};

  double result = 0.;
  double exp    = 1.;
  for (usize i = 0; i <= order; ++i)
  {
    result += coeff[i] * exp;
    exp    *= x;
  }

  return result;
}


double answers[] =
{
  +9.86605780000000000e-1,
  +9.86605780000000000e-1,
  +1.11205028666667000e+0,
  +1.11205028666667000e+0,
  +1.18415219466667000e+0,
  +1.18415219466667000e+0,
  +1.19541589752381000e+0,
  +1.19541589752381000e+0,
  +1.25352482641270000e+0,
  +1.25352482641270000e+0,
  +1.33934834641270000e+0,
  +1.33934834641270000e+0,
  +1.37858423256654000e+0,
  +1.37858423256654000e+0,
  +1.44621663923321000e+0,
  +1.44621663923321000e+0,
};


int test_quadrature()
{
  for (usize o = 0; o <= 15; ++o)
  {
    // find smallest n to integrate polynomial of order o exactly
    usize n = (usize)ceil((double(o) + 1.) / 2.);

    quadrature_rule rule = gauss_legendre_rules[n - 1];

    double integral = 0.;
    for (usize i = 0; i < rule.n; ++i)
    {
      integral += polynomial(o, rule.x[i]) * rule.w[i];
    }

    if ((answers[o] - integral) > 1e-14)
    {
      printf("quadrature test failed at order %zu!\n", o);
      exit(1);
    }
  }

  return 0;
}

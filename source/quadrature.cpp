#pragma once


#include <cstdio>
#include <cstdlib>

#include "types.cpp"


struct quadrature_rule
{
  usize         n;
  const double* x;
  const double* w;
};


const usize  gauss_legendre_n1_n                      = 1;
const double gauss_legendre_n1_x[gauss_legendre_n1_n] = {+0.00000000000000000e+00};
const double gauss_legendre_n1_w[gauss_legendre_n1_n] = {+2.00000000000000000e+00};

const usize  gauss_legendre_n2_n                      = 2;
const double gauss_legendre_n2_x[gauss_legendre_n2_n] = {-5.77350269189625731e-01, +5.77350269189625731e-01};
const double gauss_legendre_n2_w[gauss_legendre_n2_n] = {+1.00000000000000000e+00, +1.00000000000000000e+00};

const usize  gauss_legendre_n3_n                      = 3;
const double gauss_legendre_n3_x[gauss_legendre_n3_n] = {-7.74596669241483404e-01, +0.00000000000000000e+00, +7.74596669241483404e-01};
const double gauss_legendre_n3_w[gauss_legendre_n3_n] = {+5.55555555555555802e-01, +8.88888888888888284e-01, +5.55555555555555802e-01};

const usize  gauss_legendre_n4_n                      = 4;
const double gauss_legendre_n4_x[gauss_legendre_n4_n] = {-8.61136311594052573e-01, -3.39981043584856313e-01, +3.39981043584856313e-01, +8.61136311594052573e-01};
const double gauss_legendre_n4_w[gauss_legendre_n4_n] = {+3.47854845137453794e-01, +6.52145154862546206e-01, +6.52145154862546206e-01, +3.47854845137453794e-01};

const usize  gauss_legendre_n5_n                      = 5;
const double gauss_legendre_n5_x[gauss_legendre_n5_n] = {-9.06179845938663964e-01, -5.38469310105683108e-01, +0.00000000000000000e+00, +5.38469310105683108e-01, +9.06179845938663964e-01};
const double gauss_legendre_n5_w[gauss_legendre_n5_n] = {+2.36926885056188974e-01, +4.78628670499366526e-01, +5.68888888888888999e-01, +4.78628670499366526e-01, +2.36926885056188974e-01};

const usize  gauss_legendre_n6_n                      = 6;
const double gauss_legendre_n6_x[gauss_legendre_n6_n] = {-9.32469514203152050e-01, -6.61209386466264482e-01, -2.38619186083196932e-01, +2.38619186083196932e-01, +6.61209386466264482e-01, +9.32469514203152050e-01};
const double gauss_legendre_n6_w[gauss_legendre_n6_n] = {+1.71324492379170162e-01, +3.60761573048138551e-01, +4.67913934572691204e-01, +4.67913934572691204e-01, +3.60761573048138551e-01, +1.71324492379170162e-01};

const usize  gauss_legendre_n7_n                      = 7;
const double gauss_legendre_n7_x[gauss_legendre_n7_n] = {-9.49107912342758375e-01, -7.41531185599394460e-01, -4.05845151377397184e-01, +0.00000000000000000e+00, +4.05845151377397184e-01, +7.41531185599394460e-01, +9.49107912342758375e-01};
const double gauss_legendre_n7_w[gauss_legendre_n7_n] = {+1.29484966168869953e-01, +2.79705391489276700e-01, +3.81830050505118757e-01, +4.17959183673469070e-01, +3.81830050505118757e-01, +2.79705391489276700e-01, +1.29484966168869953e-01};

const usize  gauss_legendre_n8_n                      = 8;
const double gauss_legendre_n8_x[gauss_legendre_n8_n] = {-9.60289856497536287e-01, -7.96666477413626728e-01, -5.25532409916328991e-01, -1.83434642495649836e-01, +1.83434642495649836e-01, +5.25532409916328991e-01, +7.96666477413626728e-01, +9.60289856497536287e-01};
const double gauss_legendre_n8_w[gauss_legendre_n8_n] = {+1.01228536290375620e-01, +2.22381034453374732e-01, +3.13706645877887602e-01, +3.62683783378362046e-01, +3.62683783378362046e-01, +3.13706645877887602e-01, +2.22381034453374732e-01, +1.01228536290375620e-01};


const usize num_gauss_legendre_rules = 8;
quadrature_rule gauss_legendre_rules[num_gauss_legendre_rules] = {
  {gauss_legendre_n1_n, gauss_legendre_n1_x, gauss_legendre_n1_w},
  {gauss_legendre_n2_n, gauss_legendre_n2_x, gauss_legendre_n2_w},
  {gauss_legendre_n3_n, gauss_legendre_n3_x, gauss_legendre_n3_w},
  {gauss_legendre_n4_n, gauss_legendre_n4_x, gauss_legendre_n4_w},
  {gauss_legendre_n5_n, gauss_legendre_n5_x, gauss_legendre_n5_w},
  {gauss_legendre_n6_n, gauss_legendre_n6_x, gauss_legendre_n6_w},
  {gauss_legendre_n7_n, gauss_legendre_n7_x, gauss_legendre_n7_w},
  {gauss_legendre_n8_n, gauss_legendre_n8_x, gauss_legendre_n8_w},
};


// Picks Gauss-Legendre quadrature rule based on element order.
// If we're looking for quadrature rules accurate to order 2 * p + 1, and
// Gauss-Legendre quadrature rules integrate exactly for polynomials up to order
// 2 * n - 1, then the n selected is always p + 1.
//
// By default, this function uses the number of points required to integrate
// polynomials of order 2 * p + 1 exactly, but extra points may be optionally
// requested.

quadrature_rule get_quadrature_rule(usize p, usize extra_points = 0)
{
  usize n        = p + 1 + extra_points;
  usize rule_num = n - 1;

  if (rule_num > num_gauss_legendre_rules)
  {
    printf("requested quadrature rule with %zu points, \n"
           "  only quadrature rules up to %zu points supported!\n",
           n, num_gauss_legendre_rules);
    exit(1);
  }

  return gauss_legendre_rules[rule_num];
}

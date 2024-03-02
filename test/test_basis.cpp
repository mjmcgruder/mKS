#pragma once


#include <cstdio>

#include "types.cpp"
#include "basis.cpp"
#include "quadrature.cpp"


int test_basis()
{
  const usize     p_max = 7;
  const double    tol   = 1e-14;
  quadrature_rule qrule = get_quadrature_rule(p_max);

  double b[p_max + 1];

  bool fail_flag = false;
  for (usize p = 0; p < p_max; ++p)
  {
    for (usize iq = 0; iq < qrule.n; ++iq)
    {
      lagrange_chebyshev_bf(p, qrule.x[iq], b);

      double sum = 0.;
      for (usize ib = 0; ib < p + 1; ++ib)
      {
        sum += b[ib];
      }

      double err = std::abs(sum - 1.);
      if (err > tol)
      {
        printf(
        "Basis function summation gave %+.17e, should be 1! (off by %+.17e)\n",
        sum, err);
        fail_flag = true;
      }
    }
  }

  return (fail_flag) ? 1 : 0;
}


int test_basis_d1()
{
  const usize     p_max = 7;
  const double    tol   = 1e-9;  // hmmmm, derivative error increasis with order
  quadrature_rule qrule = get_quadrature_rule(p_max);

  double db[p_max + 1];

  bool fail_flag = false;
  for (usize p = 0; p < p_max; ++p)
  {
    for (usize iq = 0; iq < qrule.n; ++iq)
    {
      lagrange_chebyshev_d1(p, qrule.n, db);

      double sum = 0.;
      for (usize ib = 0; ib < p + 1; ++ib)
      {
        sum += db[ib];
      }

      double err = std::abs(sum);
      if (err > tol)
      {
        printf("Basis derivative summation gave %+.17e, should be 0! (off by "
               "%+.17e)\n",
               sum, err);
        fail_flag = true;
      }
    }
  }

  return (fail_flag) ? 1 : 0;
}

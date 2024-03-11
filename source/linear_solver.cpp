#pragma once


#include <cstring>
#include <cstdio>

#include "linalg.cpp"


// TODO: b needs to be norm 1?
void GMRES(void (*matvec)(double* v, double* mv), double* b, usize m, double* x,
           usize ninner, usize nouter, double tol)
{
  // "solution" storage
  double* Qmat  = new double[m * (ninner + 1)];
  double* Hmat  = new double[(ninner + 1) * ninner];
  double* cs    = new double[ninner];
  double* sn    = new double[ninner];
  double* betae = new double[ninner + 1];
  double* y     = new double[ninner + 1];

  // "workspace" storage
  double* Ax = new double[m];
  double* r0 = new double[m];
  double* ki = new double[m];
  double* ko = new double[m];

#define Q(r, c) (Qmat[(r) * (ninner + 1) + (c)])
#define H(r, c) (Hmat[(r) * (ninner) + (c)])

  bool isconverged = false;

  for (usize o = 0; o < nouter; ++o)
  {
    memset(Qmat, 0, m * (ninner + 1)      * sizeof(double));
    memset(Hmat, 0, (ninner + 1) * ninner * sizeof(double));
    memset(cs,   0, ninner                * sizeof(double));
    memset(sn,   0, ninner                * sizeof(double));

    // first normalized krylov vector
    matvec(x, Ax);
    for (usize i = 0; i < m; ++i) { r0[i] = b[i] - Ax[i]; }
    double r0_norm = norm(r0, m);
    for (usize i = 0; i < m; ++i) { Q(i, 0) = r0[i] / r0_norm; }

    // RHS setup
    memset(betae, 0, (ninner + 1) * sizeof(double));
    betae[0] = r0_norm;

    usize n = 0;
    for (n = 0; n < ninner; ++n)
    {
      /* Arnoldi */

      for (usize i = 0; i < m; ++i) { ki[i] = Q(i, n); }
      matvec(ki, ko);  // copies for mem ordering, hmmmmm....
      for (usize i = 0; i < m; ++i) { Q(i, n + 1) = ko[i]; }

      for (usize j = 0; j < n + 1; ++j)
      {
        H(j, n) = 0.;
        for (usize i = 0; i < m; ++i) { H(j, n) += Q(i, j) * Q(i, n + 1); }
        for (usize i = 0; i < m; ++i) { Q(i, n + 1) -= H(j, n) * Q(i, j); }
      }
      for (usize i = 0; i < m; ++i) { ko[i] = Q(i, n + 1); }
      H(n + 1, n) = norm(ko, m);  // also copying for mem ordering...
      for (usize i = 0; i < m; ++i) { Q(i, n + 1) /= H(n + 1, n); }

      /* Givens */

      // previous rotations
      for (usize pn = 0; pn < n; ++pn)
      {
        double tmp   = H(pn, n) *  cs[pn] + H(pn + 1, n) * sn[pn];
        H(pn + 1, n) = H(pn, n) * -sn[pn] + H(pn + 1, n) * cs[pn];
        H(pn, n)     = tmp;
      }

      // current rotation
      double r  = H(n, n);
      double s  = H(n + 1, n);
      double id = 1. / sqrt(r * r + s * s);
      cs[n]     = r * id;
      sn[n]     = s * id;

      H(n, n)     = H(n, n) * cs[n] + H(n + 1, n) * sn[n];
      H(n + 1, n) = 0.;

      // RHS rotation
      betae[n + 1] = -sn[n] * betae[n];
      betae[n]     =  cs[n] * betae[n];

      double residual = fabs(betae[n + 1]);

      printf("%zu %zu : %+.17e\n", o, n, residual);

      if (residual < tol)
      {
        isconverged = true;
        break;
      }
    }

    if (n > ninner - 1) { n = ninner - 1; }

    for (usize r = n; r-- > 0;)
    {
      double sum = 0.;
      for (usize c = r + 1; c < n + 1; ++c)
      {
        sum += y[c] * H(r, c);
      }
      y[r] = (-sum + betae[r]) / H(r, r);
    }

    for (usize i = 0; i < m; ++i)
    {
      for (usize k = 0; k < n + 1; ++k)
      {
        x[i] += y[k] * Q(i, k);
      }
    }

    if (isconverged) { break; }
  }

#undef Q
#undef H

  delete[] Qmat;
  delete[] Hmat;
  delete[] cs;
  delete[] sn;
  delete[] betae;
  delete[] y;

  delete[] Ax;
  delete[] r0;
  delete[] ki;
  delete[] ko;
}

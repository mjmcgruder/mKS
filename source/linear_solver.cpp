#pragma once


#include <cstring>
#include <cstdio>

#include "linalg.cpp"


struct GMRES_settings
{
  usize  ninner;
  usize  nouter;
  double rtol;
};


struct GMRES_workspace
{
  // "solution" storage
  double* Qmat;
  double* Hmat;
  double* cs;
  double* sn;
  double* betae;
  double* y;

  // "workspace" storage
  double* Ax;
  double* r0;
  double* ki;
  double* ko;

  // ---

  GMRES_workspace(GMRES_settings settings, usize m);
  ~GMRES_workspace();

  GMRES_workspace(const GMRES_workspace& oth)            = delete;
  GMRES_workspace(GMRES_workspace&& oth)                 = delete;
  GMRES_workspace& operator=(const GMRES_workspace& oth) = delete;
  GMRES_workspace& operator=(GMRES_workspace&& oth)      = delete;
};


GMRES_workspace::GMRES_workspace(GMRES_settings settings, usize m)
{
  Qmat  = new double[m * (settings.ninner + 1)];
  Hmat  = new double[(settings.ninner + 1) * settings.ninner];
  cs    = new double[settings.ninner];
  sn    = new double[settings.ninner];
  betae = new double[settings.ninner + 1];
  y     = new double[settings.ninner + 1];

  Ax = new double[m];
  r0 = new double[m];
  ki = new double[m];
  ko = new double[m];
}


GMRES_workspace::~GMRES_workspace()
{
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


// Solves Ax = b using GMRES until either the specified residual tolerance or
// one of the maximum iteration counts is met. "A" is supplied as a functor
// allowing matrix free versions of the algorithm.

template<typename T>
void GMRES(T& A, double* b, double* x, usize m,
           GMRES_settings& set, GMRES_workspace& wsp)
{
#define Q(r, c) (wsp.Qmat[(r) * (set.ninner + 1) + (c)])
#define H(r, c) (wsp.Hmat[(r) * (set.ninner) + (c)])

  bool isconverged = false;

  for (usize o = 0; o < set.nouter; ++o)
  {
    memset(wsp.Qmat, 0, m * (set.ninner + 1)          * sizeof(double));
    memset(wsp.Hmat, 0, (set.ninner + 1) * set.ninner * sizeof(double));
    memset(wsp.cs,   0, set.ninner                    * sizeof(double));
    memset(wsp.sn,   0, set.ninner                    * sizeof(double));

    // first normalized krylov vector
    A.matvec(x, wsp.Ax);
    for (usize i = 0; i < m; ++i) { wsp.r0[i] = b[i] - wsp.Ax[i]; }
    double r0_norm = norm(wsp.r0, m);
    for (usize i = 0; i < m; ++i) { Q(i, 0) = wsp.r0[i] / r0_norm; }

    // RHS setup
    memset(wsp.betae, 0, (set.ninner + 1) * sizeof(double));
    wsp.betae[0] = r0_norm;

    usize n = 0;
    for (n = 0; n < set.ninner; ++n)
    {
      /* Arnoldi */

      for (usize i = 0; i < m; ++i) { wsp.ki[i] = Q(i, n); }
      A.matvec(wsp.ki, wsp.ko);  // copies for mem ordering, hmmmmm....
      for (usize i = 0; i < m; ++i) { Q(i, n + 1) = wsp.ko[i]; }

      for (usize j = 0; j < n + 1; ++j)
      {
        H(j, n) = 0.;
        for (usize i = 0; i < m; ++i) { H(j, n) += Q(i, j) * Q(i, n + 1); }
        for (usize i = 0; i < m; ++i) { Q(i, n + 1) -= H(j, n) * Q(i, j); }
      }
      for (usize i = 0; i < m; ++i) { wsp.ko[i] = Q(i, n + 1); }
      H(n + 1, n) = norm(wsp.ko, m);  // also copying for mem ordering...
      for (usize i = 0; i < m; ++i) { Q(i, n + 1) /= H(n + 1, n); }

      /* Givens */

      // previous rotations
      for (usize pn = 0; pn < n; ++pn)
      {
        double tmp   = H(pn, n) *  wsp.cs[pn] + H(pn + 1, n) * wsp.sn[pn];
        H(pn + 1, n) = H(pn, n) * -wsp.sn[pn] + H(pn + 1, n) * wsp.cs[pn];
        H(pn, n)     = tmp;
      }

      // current rotation
      double r  = H(n, n);
      double s  = H(n + 1, n);
      double id = 1. / sqrt(r * r + s * s);
      wsp.cs[n] = r * id;
      wsp.sn[n] = s * id;

      H(n, n)     = H(n, n) * wsp.cs[n] + H(n + 1, n) * wsp.sn[n];
      H(n + 1, n) = 0.;

      // RHS rotation
      wsp.betae[n + 1] = -wsp.sn[n] * wsp.betae[n];
      wsp.betae[n]     =  wsp.cs[n] * wsp.betae[n];

      double residual = fabs(wsp.betae[n + 1]);

      // printf("  %zu %zu : %+.17e\n", o, n, residual);

      if (residual < set.rtol * r0_norm)
      {
        isconverged = true;
        break;
      }
    }

    if (n > set.ninner - 1) { n = set.ninner - 1; }

    // triangular solve
    for (usize r = n + 1; r-- > 0;)
    {
      double sum = 0.;
      for (usize c = r + 1; c < n + 1; ++c)
      {
        sum += wsp.y[c] * H(r, c);
      }
      wsp.y[r] = (-sum + wsp.betae[r]) / H(r, r);
    }

    // solution update
    for (usize i = 0; i < m; ++i)
    {
      for (usize k = 0; k < n + 1; ++k)
      {
        x[i] += wsp.y[k] * Q(i, k);
      }
    }

    A.matvec(x, wsp.ko);  // ko not used at the moment
    for (usize i = 0; i < m; ++i) { wsp.ko[i] -= b[i]; }
    double real_resid = norm(wsp.ko, m);
    printf("  real resid: %+.17e\n", real_resid);

    if (isconverged) { break; }
  }

#undef Q
#undef H
}

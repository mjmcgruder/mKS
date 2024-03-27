#pragma once


#include <cfloat>

#include "residual.cpp"
#include "geometry.cpp"
#include "types.cpp"
#include "linalg.cpp"


// Finds an appropriate perturbed state for finite difference based matrix
// vector products. Returns the step size used if successful, returns DBL_MAX if
// the norm of the perturbation direction "v" is zero.

double perturb_state(const double* uk, const double* v, const usize m, double* up)
{
  double l2vnorm = l2norm(v, m);
  if (l2vnorm < DBL_EPSILON)
  {
    return DBL_MAX;
  }

  // step size calculation from:
  // A Jacobian-free approximate Newton–Krylov startup strategy for RANS simulations

  const double erel = 1e-6;
  const double umin = 1e-4;

  double l1vnorm = l1norm(v, m);
  double vdotu   = dot(uk, v, m);

  double eps;
  if (std::abs(vdotu) > umin * l1vnorm)
  {
    eps = erel * vdotu / (l2vnorm * l2vnorm);
  }
  else
  {
    double sign = (vdotu >= 0.) ? +1. : -1.;
    eps = erel * umin * sign * l1vnorm / (l2vnorm * l2vnorm);
  }

  for (usize i = 0; i < m; ++i) { up[i] = uk[i] + eps * v[i]; }

  return eps;
}


class BDF1
{
  double*    uk;     // current nonlinear solve guess (not owned)

  // previous state (multi step method)
  double*    un;     // previous step
  // solution information
  comp_geom& geom;   // computational geometry
  // perturbed workspace
  double*    up;     // perturbed linearization state, for finite difference
  double*    Fp;     // perturbed state unsteady residual
  // steady residual worksapce
  double*    resid;  // residual space
  double*    f;      // steady residual

 public:

  double* F;   // function eval for newton and Jacobian
  double  dt;  // time step, available for pseudo transient continuation

  // computational routines

  void uresid(double* u, double* F);   // (unsteady) resiudal
  void matvec(double* v, double* mv);  // jacobian vector product

  // internal state update routines

  void update_guess(double* st);
  void update_step(double* st);

  // --

  BDF1(comp_geom& geom, double dt);
  ~BDF1();

  BDF1(const BDF1& oth)            = delete;
  BDF1(BDF1&& oth)                 = delete;
  BDF1& operator=(const BDF1& oth) = delete;
  BDF1& operator=(BDF1&& oth)      = delete;
};


void BDF1::uresid(double* u, double* F)
{
  compute_residual(geom, u, resid, f);

  for (usize i = 0; i < geom.size_st(); ++i)
  {
    F[i] = u[i] - un[i] - dt * f[i];
  }
}


void BDF1::matvec(double* v, double* mv)
{
  double eps = perturb_state(uk, v, geom.size_st(), up);
  if (eps == DBL_MAX)
  {
    memset(mv, 0, geom.size_st() * sizeof(double));
    return;
  }

  uresid(up, Fp);

  double ieps = 1. / eps;
  for (usize i = 0; i < geom.size_st(); ++i) { mv[i] = (Fp[i] - F[i]) * ieps; }
}


void BDF1::update_guess(double* new_uk)
{
  uk = new_uk;
  uresid(uk, F);
}


void BDF1::update_step(double* u)
{
  memcpy(un, u, geom.size_st() * sizeof(*u));
}


BDF1::BDF1(comp_geom& geom_, double dt_) : geom(geom_), dt(dt_)
{
  un    = new double[geom.size_st()];
  resid = new double[geom.size_st()];
  f     = new double[geom.size_st()];
  up    = new double[geom.size_st()];
  F     = new double[geom.size_st()];
  Fp    = new double[geom.size_st()];
}


BDF1::~BDF1()
{
  delete[] un;
  delete[] resid;
  delete[] f;
  delete[] up;
  delete[] F;
  delete[] Fp;
}


class BDF2
{
  double*    uk;     // current nonlinear solve guess (not owned)

  // previous state (multi step method)
  double*    un;     // previous step
  double*    unm1;   // previous previous step
  // solution information
  comp_geom& geom;   // computational geometry
  // perturbed workspace
  double*    up;     // perturbed linearization state, for finite difference
  double*    Fp;     // perturbed state unsteady residual
  // steady residual worksapce
  double*    resid;  // residual space
  double*    f;      // steady residual

 public:

  double* F;   // function eval for newton and Jacobian
  double  dt;  // time step, available for pseudo transient continuation

  // computational routines

  void uresid(double* u, double* F);   // (unsteady) resiudal
  void matvec(double* v, double* mv);  // jacobian vector product

  // internal state update routines

  void update_guess(double* st);
  void update_step(double* st);

  // --

  BDF2(comp_geom& geom, double dt);
  ~BDF2();

  BDF2(const BDF2& oth)            = delete;
  BDF2(BDF2&& oth)                 = delete;
  BDF2& operator=(const BDF2& oth) = delete;
  BDF2& operator=(BDF2&& oth)      = delete;
};


void BDF2::uresid(double* u, double* F)
{
  compute_residual(geom, u, resid, f);

  for (usize i = 0; i < geom.size_st(); ++i)
  {
    F[i] = 0.5 * unm1[i] - 2. * un[i] + 1.5 * u[i] - dt * f[i];
  }
}


void BDF2::matvec(double* v, double* mv)
{
  double eps = perturb_state(uk, v, geom.size_st(), up);
  if (eps == DBL_MAX)
  {
    memset(mv, 0, geom.size_st() * sizeof(double));
    return;
  }

  uresid(up, Fp);

  double ieps = 1. / eps;
  for (usize i = 0; i < geom.size_st(); ++i) { mv[i] = (Fp[i] - F[i]) * ieps; }
}


void BDF2::update_guess(double* new_uk)
{
  uk = new_uk;
  uresid(uk, F);
}


void BDF2::update_step(double* u)
{
  memcpy(unm1, un, geom.size_st() * sizeof(*un));
  memcpy(un, u, geom.size_st() * sizeof(*u));
}


BDF2::BDF2(comp_geom& geom_, double dt_) : geom(geom_), dt(dt_)
{
  un    = new double[geom.size_st()];
  unm1  = new double[geom.size_st()];
  resid = new double[geom.size_st()];
  f     = new double[geom.size_st()];
  up    = new double[geom.size_st()];
  F     = new double[geom.size_st()];
  Fp    = new double[geom.size_st()];
}


BDF2::~BDF2()
{
  delete[] un;
  delete[] unm1;
  delete[] resid;
  delete[] f;
  delete[] up;
  delete[] F;
  delete[] Fp;
}

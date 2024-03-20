#pragma once


#include <cfloat>

#include "residual.cpp"
#include "geometry.cpp"
#include "types.cpp"
#include "linalg.cpp"


class BDF1
{
  double*    uk;     // current nonlinear solve guess (not owned)

  // previous state (multi step method)
  double*    un;     // previous step
  // solution information
  comp_geom& geom;   // computational geometry
  double     dt;     // time step
  // perturbed workspace
  double*    up;     // perturbed linearization state, for finite difference
  double*    Fp;     // perturbed state unsteady residual
  // steady residual worksapce
  double*    resid;  // residual space
  double*    f;      // steady residual

 public:

  double* F;  // function eval for newton and Jacobian

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
  const double erel = 1e-4;

  double vnorm = norm(v, geom.size_st());
  if (vnorm < DBL_EPSILON)
  {
    memset(mv, 0, geom.size_st() * sizeof(double));
    return;
  }

  double vdotu = std::abs(dot(uk, v, geom.size_st()));

  double eps  = (erel / vnorm) * vdotu;
  double ieps = 1. / eps;

  for (usize i = 0; i < geom.size_st(); ++i) { up[i] = uk[i] + eps * v[i]; }

  uresid(up, Fp);

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

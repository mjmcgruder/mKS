#pragma once


#include <cfloat>

#include "residual.cpp"
#include "geometry.cpp"
#include "types.cpp"
#include "linalg.cpp"


template<typename T>
struct matvec_fd
{
  T& system;   // reference to object providing uresid routine

  usize   m;   // vector size

  double* uk;  // current linearization state (not owned)
  double* F;   // current function eval (at uk)

  double* up;  // perturbed linearization state, for finite difference
  double* Fp;  // perturbed state unsteady residual

  void matvec(double* v, double* mv);
  void update_guess(double* st);

  matvec_fd(T& system_, usize m_);
  ~matvec_fd();

  matvec_fd(const matvec_fd& oth)            = delete;
  matvec_fd(matvec_fd&& oth)                 = delete;
  matvec_fd& operator=(const matvec_fd& oth) = delete;
  matvec_fd& operator=(matvec_fd&& oth)      = delete;
};


template<typename T>
void matvec_fd<T>::matvec(double* v, double* mv)
{
  double l2vnorm = l2norm(v, m);
  if (l2vnorm < DBL_EPSILON)
  {
    memset(mv, 0, m * sizeof(double));
    return;
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

  system.uresid(up, Fp);

  double ieps = 1. / eps;
  for (usize i = 0; i < m; ++i) { mv[i] = (Fp[i] - F[i]) * ieps; }
}


template<typename T>
void matvec_fd<T>::update_guess(double* new_uk)
{
  uk = new_uk;
  system.uresid(uk, F);
}


template<typename T>
matvec_fd<T>::matvec_fd(T& system_, usize m_) : system(system_), m(m_)
{
  F  = new double[m];
  up = new double[m];
  Fp = new double[m];
}


template<typename T>
matvec_fd<T>::~matvec_fd()
{
  delete[] F;
  delete[] up;
  delete[] Fp;
}


template<usize nstep>
struct multi_step_store
{
  comp_geom& geom;

  double* unmx[nstep];  // stores most recent previous step last

  double* resid;
  double* f;

  void update(double* u)
  {
    for (usize i = 0; i < nstep - 1; ++i)
    {
      memcpy(unmx[i], unmx[i + 1], geom.size_st() * sizeof(*u));
    }
    memcpy(unmx[nstep - 1], u, geom.size_st() * sizeof(*u));
  }

  multi_step_store(comp_geom& geom_) : geom(geom_)
  {
    for (usize i = 0; i < nstep; ++i)
    {
      unmx[i] = new double[geom.size_st()];
    }

    resid = new double[geom.size_st()];
    f     = new double[geom.size_st()];
  }

  ~multi_step_store()
  {
    for (usize i = 0; i < nstep; ++i)
    {
      delete[] unmx[i];
    }

    delete[] resid;
    delete[] f;
  }

  multi_step_store(const multi_step_store& oth)            = delete;
  multi_step_store(multi_step_store&& oth)                 = delete;
  multi_step_store& operator=(const multi_step_store& oth) = delete;
  multi_step_store& operator=(multi_step_store&& oth)      = delete;
};


struct time_scheme_steady
{
  comp_geom& geom;

  double* resid;

  void uresid(double* u, double* F)
  {
    compute_residual(geom, u, resid, F);
  }

  time_scheme_steady(comp_geom& geom_) : geom(geom_)
  {
    resid = new double[geom.size_st()];
  }
  ~time_scheme_steady()
  {
    delete[] resid;
  }
};


struct BDF1
{
  multi_step_store<1> store;
  double  dt;  // time step, available for pseudo transient continuation

  void uresid(double* u, double* F)
  {
    compute_residual(store.geom, u, store.resid, store.f);

    for (usize i = 0; i < store.geom.size_st(); ++i)
    {
      F[i] = u[i] - store.unmx[0][i] - dt * store.f[i];
    }
  }

  void update_step(double* u)
  {
    store.update(u);
  }

  BDF1(comp_geom& geom_, double dt_) : store(geom_), dt(dt_)
  {}
};


struct BDF2
{
  multi_step_store<2> store;
  double  dt;  // time step, available for pseudo transient continuation

  void uresid(double* u, double* F)
  {
    compute_residual(store.geom, u, store.resid, store.f);

    for (usize i = 0; i < store.geom.size_st(); ++i)
    {
      F[i] = 0.5 * store.unmx[0][i] - 2. * store.unmx[1][i] + 1.5 * u[i] - dt * store.f[i];
    }
  }

  void update_step(double* u)
  {
    store.update(u);
  }

  BDF2(comp_geom& geom_, double dt_) : store(geom_), dt(dt_)
  {}
};


const double alpha = (2. / sqrt(3.)) * cos(PI / 18.);


struct DIRK_wrksp
{
  const double* yn;  // (not owned) pointer to current state

  comp_geom&    geom;
  double        dt;

  double*       u;
  double*       resid;
  double*       f;

  DIRK_wrksp(comp_geom& geom_, double dt_) : geom(geom_), dt(dt_)
  {
    u     = new double[geom.size_st()];
    resid = new double[geom.size_st()];
    f     = new double[geom.size_st()];
  }

  ~DIRK_wrksp()
  {
    delete[] u;
    delete[] resid;
    delete[] f;
  }
};


struct DIRK_1
{
  DIRK_wrksp& wrk;

  void uresid(double* k1, double* F)
  {
    for (usize i = 0; i < wrk.geom.size_st(); ++i)
    {
      wrk.u[i] = wrk.yn[i] + wrk.dt * 0.5 * (1. + alpha) * k1[i];
    }

    compute_residual(wrk.geom, wrk.u, wrk.resid, wrk.f);

    for (usize i = 0; i < wrk.geom.size_st(); ++i)
    {
      F[i] = k1[i] - wrk.f[i];
    }
  }

  DIRK_1(DIRK_wrksp& wrk_) : wrk(wrk_)
  {}
};

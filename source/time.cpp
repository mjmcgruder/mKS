#pragma once


#include <cmath>
#include <cfloat>

#include "linalg.cpp"
#include "residual.cpp"
#include "linear_solver.cpp"


void euler(comp_geom& geom, usize nt, double dt, double* state)
{
  double* resid = new double[geom.size_st()];
  double* f     = new double[geom.size_st()];

  for (usize it = 0; it < nt; ++it)
  {
    compute_residual(geom, state, resid, f);

    printf("%zu: %+.3e\n", it, norm_perelem(f, geom.size_st()));

    for (usize i = 0; i < geom.size_st(); ++i)
    {
      state[i] += dt * f[i];
    }
  }

  delete[] resid;
  delete[] f;
}


void rk4(comp_geom& geom, usize nt, double dt, double* y1, double* hist)
{
  double* resid = new double[geom.size_st()];

  double* k1 = new double[geom.size_st()];
  double* k2 = new double[geom.size_st()];
  double* k3 = new double[geom.size_st()];
  double* k4 = new double[geom.size_st()];

  double* y2 = new double[geom.size_st()];
  double* y3 = new double[geom.size_st()];
  double* y4 = new double[geom.size_st()];

  for (usize it = 0; it < nt; ++it)
  {
    memcpy(hist + it * geom.size_st(), y1, geom.size_st() * sizeof(double));

    compute_residual(geom, y1, resid, k1);

    printf("%zu: %+.3e\n", it, norm_perelem(k1, geom.size_st()));

    for (usize i = 0; i < geom.size_st(); ++i) y2[i] = y1[i] + 0.5 * dt * k1[i];

    compute_residual(geom, y2, resid, k2);

    for (usize i = 0; i < geom.size_st(); ++i) y3[i] = y1[i] + 0.5 * dt * k2[i];

    compute_residual(geom, y3, resid, k3);

    for (usize i = 0; i < geom.size_st(); ++i) y4[i] = y1[i] + dt * k3[i];

    compute_residual(geom, y4, resid, k4);

    for (usize i = 0; i < geom.size_st(); ++i)
    {
      y1[i] += (dt / 6.) * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    }
  }

  memcpy(hist + nt * geom.size_st(), y1, geom.size_st() * sizeof(double));

  delete[] resid;

  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] k4;

  delete[] y2;
  delete[] y3;
  delete[] y4;
}


// // TODO: in future, template this and provide time scheme specific
// // specializations, maybe a value template with time scheme enum or something
// struct BDF1
// {
//   // simulation information
//   comp_geom& geom;
//   double     dt;
//
//   // known time steps (owned)
//   double* un;
//
//   // linearization state (current nonlinear solve guess) (not owned)
//   double* uk;
//
//   // residual work space (owned)
//   double* resid;
//   double* f;        // steady residual
//   double* uk_pert;  // perturbed linearization state, for finite difference
//
//   double* Fo;       // original (linearization state) unsteady residual
//   double* Fp;       // perturbed state unsteady residual
//
//   // computational routines
//
//   void uresid(double* u, double* F);   // (unsteady) resiudal
//   void matvec(double* v, double* mv);  // jacobian vector product
//
//   // internal state update routines
//
//   void update_linearization_state(double* st);
//   void update_time_step(double* st);
//
//   // --
//
//   BDF1(comp_geom& geom, double dt);
//   ~BDF1();
//
//   BDF1(const BDF1& oth)            = delete;
//   BDF1(BDF1&& oth)                 = delete;
//   BDF1& operator=(const BDF1& oth) = delete;
//   BDF1& operator=(BDF1&& oth)      = delete;
// };
//
//
// void BDF1::uresid(double* u, double* F)
// {
//   compute_residual(geom, u, resid, f);
//
//   for (usize i = 0; i < geom.size_st(); ++i)
//   {
//     F[i] = u[i] - un[i] - dt * f[i];
//   }
// }
//
//
// void BDF1::matvec(double* v, double* mv)
// {
//   double normv = norm(v, geom.size_st()) + DBL_EPSILON;
//   double eps   = 1e-4;
//
//   for (usize i = 0; i < geom.size_st(); ++i)
//   {
//     uk_pert[i] = uk[i] + eps * (v[i] / normv);
//   }
//
//   uresid(uk_pert, Fp);
//
//   for (usize i = 0; i < geom.size_st(); ++i)
//   {
//     mv[i] = (Fp[i] - Fo[i]) / eps;
//   }
// }
//
//
// void BDF1::update_linearization_state(double* new_uk)
// {
//   uk = new_uk;
//   uresid(uk, Fo);
// }
//
//
// void BDF1::update_time_step(double* u)
// {
//   memcpy(un, u, geom.size_st() * sizeof(*u));
// }
//
//
// BDF1::BDF1(comp_geom& geom_, double dt_) : geom(geom_), dt(dt_)
// {
//   un      = new double[geom.size_st()];
//   resid   = new double[geom.size_st()];
//   f       = new double[geom.size_st()];
//   uk_pert = new double[geom.size_st()];
//   Fo      = new double[geom.size_st()];
//   Fp      = new double[geom.size_st()];
// }
//
//
// BDF1::~BDF1()
// {
//   delete[] un;
//   delete[] resid;
//   delete[] f;
//   delete[] uk_pert;
//   delete[] Fo;
//   delete[] Fp;
// }


// void line_search_backtrack(BDF1& bdf1, double* tmpu, double* x, const double* uk)
// {
//   double norm_Fo = norm(bdf1.Fo, bdf1.geom.size_st());
//   double norm_Fp;
//   double s = 1.;
//
//   do
//   {
//     for (usize i = 0; i < bdf1.geom.size_st(); ++i) { tmpu[i] = uk[i] + x[i]; }
//
//     bdf1.uresid(tmpu, bdf1.Fp);
//     double norm_Fp = norm(bdf1.Fp, bdf1.geom.size_st());
//
//     printf("line searching: initial: %f final: %f s=%f\n", norm_Fo, norm_Fp, s);
//
//     for (usize i = 0; i < bdf1.geom.size_st(); ++i) { x[i] *= 0.5; }
//     s *= 0.5;
//   }
//   while (norm_Fp > norm_Fo);
//
//   // TODO: assign temp step to final step
// }


struct newton_workspace
{
  double* x_upd;
  double* b;

  // ---

  newton_workspace(usize m);
  ~newton_workspace();

  newton_workspace(const newton_workspace& oth)            = delete;
  newton_workspace(newton_workspace&& oth)                 = delete;
  newton_workspace& operator=(const newton_workspace& oth) = delete;
  newton_workspace& operator=(newton_workspace&& oth)      = delete;
};


newton_workspace::newton_workspace(usize m)
{
  x_upd = new double[m];
  b     = new double[m];
}


newton_workspace::~newton_workspace()
{
  delete[] x_upd;
  delete[] b;
}


struct newton_settings
{
  usize ninner;
  usize nouter;
  usize tol_lin;
  usize tol_nonlin;
};


template<typename T>
void newton_step(T& system, double* x, usize m, newton_settings& settings,
                 newton_workspace& wsp)
{
  system.update_guess(x);

  double nonlinear_residual = norm(system.Fo, m);
  printf("nonlinear residual: %+.17e\n", nonlinear_residual);

  while (nonlinear_residual > settings.tol_nonlin)
  {
    for (usize i = 0; i < m; ++i) { wsp.b[i] = -system.Fo[i]; }
    memset(wsp.x_upd, 0, m * sizeof(double));  // guess old state

    GMRES(system, wsp.b, wsp.x_upd, m,
          settings.ninner, settings.nouter, settings.tol_lin);

    for (usize i = 0; i < m; ++i) { x[i] += wsp.x_upd[i]; }  // full newton step

    system.update_guess(x);

    nonlinear_residual = norm(system.Fo, m);
    printf("nonlinear residual: %+.17e\n", nonlinear_residual);
  }
}

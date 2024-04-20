#pragma once


#include <cmath>

#include "linalg.cpp"
#include "residual.cpp"
#include "linear_solver.cpp"
#include "time_schemes.cpp"


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
  GMRES_settings gmres;
  double         nonlin_tol;
  usize          max_iter;
};


template<typename T>
void newton_solve(T& system, double* x, usize m, newton_settings& settings,
                  GMRES_workspace& gwsp, newton_workspace& nwsp)
{
  system.update_guess(x);

  double nonlinear_residual = l2norm(system.F, m);
  vrbprintf(2, "nonlinear residual: %+.8e\n", nonlinear_residual);

  usize newton_iters = 0;
  while (nonlinear_residual > settings.nonlin_tol &&
         newton_iters < settings.max_iter)
  {
    for (usize i = 0; i < m; ++i) { nwsp.b[i] = -system.F[i]; }
    memset(nwsp.x_upd, 0, m * sizeof(double));  // guess old state

    GMRES(system, nwsp.b, nwsp.x_upd, m, settings.gmres, gwsp);

    for (usize i = 0; i < m; ++i) { x[i] += nwsp.x_upd[i]; }  // full newton step

    system.update_guess(x);

    nonlinear_residual = l2norm(system.F, m);
    vrbprintf(2, "finished nonlinear iter: %zu --- residual: %+.8e\n",
              newton_iters, nonlinear_residual);

    ++newton_iters;
  }
}


enum struct evolution_scheme
{
  rk4,
  bdf2,
  count
};

const char* evolution_scheme_map[(int)evolution_scheme::count] =
{
  "rk4",
  "bdf2"
};


void evolve_rk4(comp_geom& geom, usize nt, double dt, double* y1, double* hist)
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

    vrbprintf(1, "%zu: %+.3e\n", it, l2norm(k1, geom.size_st()));

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


void evolve_bdf2(comp_geom& geom, usize nt, double dt, newton_settings& settings,
                 double* state, double* hist)
{
  newton_workspace nwsp(geom.size_st());
  GMRES_workspace  gwsp(settings.gmres, geom.size_st());

  BDF1 bdf1(geom, dt);
  BDF2 bdf2(geom, dt);

  matvec_fd<BDF1> bdf1_linsys(bdf1, geom.size_st());
  matvec_fd<BDF2> bdf2_linsys(bdf2, geom.size_st());

  bdf1.update_step(state);
  bdf2.update_step(state);

  for (usize it = 0; it < nt; ++it)
  {
    vrbprintf(2, "\n");
    vrbprintf(1, "TIME STEP %zu\n", it);
    vrbprintf(2, "\n");

    memcpy(hist + it * geom.size_st(), state, geom.size_st() * sizeof(double));

    if (it == 0)
    {
      newton_solve(bdf1_linsys, state, geom.size_st(), settings, gwsp, nwsp);
    }
    else
    {
      newton_solve(bdf2_linsys, state, geom.size_st(), settings, gwsp, nwsp);
    }

    bdf2.update_step(state);
  }

  memcpy(hist + nt * geom.size_st(), state, geom.size_st() * sizeof(double));
}

#pragma once


#include "time_schemes.cpp"
#include "time.cpp"
#include "initialization.cpp"
#include "io.cpp"


constexpr double GAUSSIAN_A = 1.;
constexpr double GAUSSIAN_B = PI;
constexpr double GAUSSIAN_C = 0.5;


double gaussian_func(double x)
{
  const double a = GAUSSIAN_A;
  const double b = GAUSSIAN_B;
  const double c = GAUSSIAN_C;

  return a * exp(-((x - b) * (x - b)) / (2. * c * c));
}


double gaussian_deriv(double x)
{
  const double a = GAUSSIAN_A;
  const double b = GAUSSIAN_B;
  const double c = GAUSSIAN_C;

  return -(a * (x - b) * exp(-((x - b) * (x - b)) / (2. * c * c))) / (c * c);
}


double manufactured_source_2(double x)
{
  const double a = GAUSSIAN_A;
  const double b = GAUSSIAN_B;
  const double c = GAUSSIAN_C;

  return (a * exp(-((x - b) * (x - b)) / (2. * c * c)) *
         (b * b - 2. * b * x - c * c + x * x)) / (c * c * c * c);
}


double manufactured_source_4(double x)
{
  const double a = GAUSSIAN_A;
  const double b = GAUSSIAN_B;
  const double c = GAUSSIAN_C;

  return -(a * exp(-((x - b) * (x - b)) / (2. * c * c)) *
         (b * b * b * b - 4. * b * b * b * x - 6. * b * b * (c * c - x * x) +
          4. * b * (3. * c * c * x - x * x * x) + 3. * c * c * c * c -
          6. * c * c * x * x + x * x * x * x)) /
         (c * c * c * c * c * c * c * c);
}


struct manufactured_result
{
  double unsteady_residual;
  double solution_error;
};


manufactured_result test_manufactured(comp_geom& geom, usize ne, usize p)
{
  manufactured_result result;

  // options

  newton_settings settings;
  settings.gmres.ninner = ne * (p + 1);
  settings.gmres.nouter = 10;
  settings.gmres.rtol   = 1e-4;
  settings.gmres.ctol   = 1e-12;
  settings.nonlin_tol   = 1e-10;
  settings.max_iter     = 20;

  // initialization

  double* state = new double[geom.size_st()];

  init_constant(geom, 0.5, state);

  time_scheme_steady steady(geom);
  matvec_fd<time_scheme_steady> steady_linsys(steady, geom.size_st());

  newton_workspace nwsp(geom.size_st());
  GMRES_workspace  gwsp(settings.gmres, geom.size_st());

  // converging

  newton_solve(steady_linsys, state, geom.size_st(), settings, gwsp, nwsp);

  // find final residual

  double* F = new double[geom.size_st()];

  // bdf1.uresid(state, F);
  steady.uresid(state, F);
  result.unsteady_residual = l2norm(F, geom.size_st());

  delete[] F;

  // calculate error ---

  quadrature_rule qrule = get_quadrature_rule(0, 15);

  double* b_eval = new double[geom.p + 1];

  double err = 0.;
  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double J   = geom.J;

    for (usize qi = 0; qi < qrule.n; ++qi)
    {
      lagrange_chebyshev_bf(geom.p, qrule.x[qi], b_eval);

      double s = 0.;
      for (usize bi = 0; bi < pp1; ++bi)
      {
        s += geom.idx_st(state, ei, bi) * b_eval[bi];
      }

      double x = geom.ref2glo(ei, qrule.x[qi]);
      double f = gaussian_func(x);

      double diff = f - s;
      double diffsq = diff * diff;

      err += diffsq * J * qrule.w[qi];
    }
  }

  result.solution_error = sqrt(err);

  delete[] b_eval;

  return result;
}


int test_convergence()
{
  double l = 2. * PI;

  printf("second order diffusion test:\n\n");

  for (usize p = 1; p <= 4; ++p)
  {
    for (usize ne = 8; ne <= 256; ne *= 2)
    {
      comp_geom geom(l, ne, p,
                     0., -1., 0., 8., 729., 9.,
                     dirichlet,
                     gaussian_func(0.), gaussian_func(l),
                     gaussian_deriv(0.), gaussian_deriv(l),
                     manufactured_source_2);

      manufactured_result res = test_manufactured(geom, ne, p);
      printf("p: %lu, ne: %3lu, resid: %.8e, error: %.8e\n",
             p, ne, res.unsteady_residual, res.solution_error);
    }
    printf("\n");
  }

  printf("fourth order diffusion test:\n\n");

  for (usize p = 2; p <= 4; ++p)
  {
    for (usize ne = 8; ne <= 128; ne *= 2)
    {
      comp_geom geom(l, ne, p,
                     0., 0., 1., 8., 729., 9.,
                     dirichlet,
                     gaussian_func(0.), gaussian_func(l),
                     gaussian_deriv(0.), gaussian_deriv(l),
                     manufactured_source_4);

      manufactured_result res = test_manufactured(geom, ne, p);
      printf("p: %lu, ne: %3lu, resid: %.8e, error: %.8e\n",
             p, ne, res.unsteady_residual, res.solution_error);
    }
    printf("\n");
  }

  printf("\n");

  return 0;
}

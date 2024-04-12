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


int test_convergence()
{
  // usize ne = 16;
  // usize p  = 2;
  //
  // double nc  = 0.;
  // double nu2 = 0.;
  // double nu4 = 1.;
  // double eta = 8.;
  // double sig = 729.;
  // double tau = 9.;
  //
  // newton_settings settings;
  // settings.gmres.ninner = ne * (p + 1);
  // settings.gmres.nouter = 10;
  // settings.gmres.rtol   = 1e-4;
  // settings.gmres.ctol   = 1e-11;
  // settings.nonlin_tol   = 1e-10;
  // settings.max_iter     = 10;
  //
  // double l = 2. * PI;
  //
  // comp_geom geom(l, ne, p,
  //                nc, nu2, nu4, eta, sig, tau,
  //                dirichlet, gaussian_func(0.), gaussian_func(l),
  //                gaussian_deriv(0.), gaussian_deriv(l),
  //                manufactured_source_4);
  //
  // double* state = new double[geom.size_st()];
  //
  // init_constant(geom, 0.5, state);
  //
  // BDF1 bdf1(geom, 0.01);
  //
  // newton_workspace nwsp(geom.size_st());
  // GMRES_workspace  gwsp(settings.gmres, geom.size_st());
  //
  // bdf1.update_step(state);
  //
  // double t = 0;
  // for (usize it = 0; it < 128; ++it)
  // {
  //   printf("\n\n");
  //   printf("TIME STEP %zu, t = %f, dt = %f\n", it, t, bdf1.dt);
  //   printf("\n\n");
  //
  //   newton_solve(bdf1, state, geom.size_st(), settings, gwsp, nwsp);
  //
  //   bdf1.update_step(state);
  //
  //   t       += bdf1.dt;
  //   if (bdf1.dt < 1.)
  //   {
  //     bdf1.dt *= 1.1;
  //   }
  // }
  //
  // // calculate error ---
  //
  // quadrature_rule qrule = get_quadrature_rule(0, 15);
  //
  // double* b_eval = new double[geom.p + 1];
  //
  // double err = 0.;
  // for (usize ei = 0; ei < geom.ne; ++ei)
  // {
  //   usize  pp1 = geom.p + 1;
  //   double J   = geom.J;
  //
  //   for (usize qi = 0; qi < qrule.n; ++qi)
  //   {
  //     lagrange_chebyshev_bf(geom.p, qrule.x[qi], b_eval);
  //
  //     double s = 0.;
  //     for (usize bi = 0; bi < pp1; ++bi)
  //     {
  //       s += geom.idx_st(state, ei, bi) * b_eval[bi];
  //     }
  //
  //     double x = geom.ref2glo(ei, qrule.x[qi]);
  //     double f = gaussian_func(x);
  //
  //     double diff = f - s;
  //     double diffsq = diff * diff;
  //
  //     err += diffsq * J * qrule.w[qi];
  //   }
  // }
  //
  // err = sqrt(err);
  //
  // delete[] b_eval;
  //
  // printf("L2 error: %.17e\n", err);
  //
  // write_hist(geom, state, 1., 1., "convergence_test.ks");
  //
  // delete[] state;

  return 0;
}

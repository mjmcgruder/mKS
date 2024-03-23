#include "initialization.cpp"
#include "io.cpp"
#include "time.cpp"
#include "time_schemes.cpp"


int main(int argc, char** argv)
{
  if (argc != 2)
  {
    printf("Usage: mKS <-h | ifile>\n");
    exit(0);
  }

  usize  ne = 8;
  usize  p  = 3;

  usize  nt = 1;
  double dt = 1e-3;

  double nc  = 1.;
  double nu2 = 1.;
  double nu4 = 0.1111;
  double eta = 16.;
  double sig = 4096.;
  double tau = 16.;

  newton_settings settings;
  settings.gmres.ninner = 100;
  settings.gmres.nouter = 5;
  settings.gmres.rtol   = 1e-3;
  settings.nonlin_tol   = 1e-8;
  settings.max_iter     = 10;

  std::map<const char*, option, option_comp> options;
  options.emplace("nelem",            mkopt(ne, "element count"));
  options.emplace("p",                mkopt(p, "polynomial order"));
  options.emplace("nstep",            mkopt(nt, "number of time steps"));
  options.emplace("dt",               mkopt(dt, "time step size"));
  options.emplace("gmres_ninner",     mkopt(settings.gmres.ninner, "max GMRES Krylov subspace size"));
  options.emplace("gmres_nouter",     mkopt(settings.gmres.nouter, "max GMRES Krylov subspace rebuilds"));
  options.emplace("gmres_rtol",       mkopt(settings.gmres.rtol, "relative linear residual convergence tolerance"));
  options.emplace("nonlin_tol",       mkopt(settings.nonlin_tol, "nonlinear residual tolerance"));
  options.emplace("nonlin_max_iters", mkopt(settings.max_iter, "max noninear iterations per time step"));
  options.emplace("eqn_nc",           mkopt(nc,  "nonlinear coefficient"));
  options.emplace("eqn_nu2",          mkopt(nu2, "second order viscosity coefficient"));
  options.emplace("eqn_nu4",          mkopt(nu4, "fourth order viscosity coefficient"));
  options.emplace("ip_eta",           mkopt(eta, "second order IP stabilization constant"));
  options.emplace("ip_sig",           mkopt(sig, "fourth order IP stabilization constant 1"));
  options.emplace("ip_tau",           mkopt(tau, "fourth order IP stabilization constant 2"));

  if (strcmp(argv[1], "-h") == 0)
  {
    for (const auto& i : options)
    {
      printf("%-20s : %s\n", i.first, i.second.desc);
    }
    exit(0);
  }
  else
  {
    read_input("ifile", options);
  }

  comp_geom geom(2. * PI, ne, p,
                 nc, nu2, nu4, eta, sig, tau,
                 periodic, 0., 0.);

  double* state = new double[geom.size_st()];
  double* hist  = new double[(nt + 1) * geom.size_st()];

  // init_gaussian(geom, 0.1, 0.5, 0.1, state);
  init_sines(geom, state);
  // init_constant(geom, 0.1, state);

  /* TIME STEPPING */

  BDF1 bdf1(geom, dt);

  newton_workspace nwsp(geom.size_st());
  GMRES_workspace  gwsp(settings.gmres, geom.size_st());

  bdf1.update_step(state);

  for (usize it = 0; it < nt; ++it)
  {
    printf("\n\n");
    printf("TIME STEP %zu\n", it);
    printf("\n\n");

    memcpy(hist + it * geom.size_st(), state, geom.size_st() * sizeof(double));

    newton_solve(bdf1, state, geom.size_st(), settings, gwsp, nwsp);

    bdf1.update_step(state);
  }

  memcpy(hist + nt * geom.size_st(), state, geom.size_st() * sizeof(double));

  write_hist(geom, hist, (nt + 1), dt, "hist.ks");

  delete[] state;
  delete[] hist;

  return 0;
}

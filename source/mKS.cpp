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

  // goemetry
  double len = 2. * PI;
  usize  ne  = 8;
  usize  p   = 3;
  // time
  usize  nt = 100;
  double dt = 1e-1;
  // equation
  double nc  = 1.;
  double nu2 = 1.;
  double nu4 = 0.1111;
  double eta = 16.;
  double sig = 4096.;
  double tau = 16.;
  // solver
  newton_settings settings;
  settings.gmres.ninner = 100;
  settings.gmres.nouter = 10;
  settings.gmres.rtol   = 1e-4;
  settings.gmres.ctol   = 1e-8;
  settings.nonlin_tol   = 1e-8;
  settings.max_iter     = 10;

  std::map<const char*, option, option_comp> options;
  options.emplace("geom_nelem",         mkopt(ne,                    "element count"));
  options.emplace("geom_p",             mkopt(p,                     "polynomial order"));
  options.emplace("time_nstep",         mkopt(nt,                    "number of time steps"));
  options.emplace("time_dt",            mkopt(dt,                    "time step size"));
  options.emplace("eqn_nc",             mkopt(nc,                    "nonlinear coefficient"));
  options.emplace("eqn_nu2",            mkopt(nu2,                   "second order viscosity coefficient"));
  options.emplace("eqn_nu4",            mkopt(nu4,                   "fourth order viscosity coefficient"));
  options.emplace("eqn_stab_eta",       mkopt(eta,                   "second order IP stabilization constant"));
  options.emplace("eqn_stab_sig",       mkopt(sig,                   "fourth order IP stabilization constant 1"));
  options.emplace("eqn_stab_tau",       mkopt(tau,                   "fourth order IP stabilization constant 2"));
  options.emplace("solve_gmres_ninner", mkopt(settings.gmres.ninner, "max GMRES Krylov subspace size"));
  options.emplace("solve_gmres_nouter", mkopt(settings.gmres.nouter, "max GMRES Krylov subspace rebuilds"));
  options.emplace("solve_gmres_rtol",   mkopt(settings.gmres.rtol,   "relative linear residual convergence tolerance"));
  options.emplace("solve_gmres_ctol",   mkopt(settings.gmres.ctol,   "cutoff linear residual tolerance"));
  options.emplace("solve_tol",          mkopt(settings.nonlin_tol,   "nonlinear residual tolerance"));
  options.emplace("solve_iters",        mkopt(settings.max_iter,     "max noninear iterations per time step"));
  options.emplace("verbosity",          mkopt(VERBOSITY,             "verbosity level [0, 3]"));

  if (strcmp(argv[1], "-h") == 0)
  {
    for (const auto& i : options)
    {
      printf("%-20s : %s\n", i.first, i.second.desc);
    }
    exit(0);
  }

  read_input(argv[1], options);

  comp_geom geom(len, ne, p,
                 nc, nu2, nu4, eta, sig, tau);

  double* state = new double[geom.size_st()];
  double* hist  = new double[(nt + 1) * geom.size_st()];

  init_sines(geom, state);

  /* TIME STEPPING */

  BDF1 bdf1(geom, dt);
  BDF2 bdf2(geom, dt);

  matvec_fd<BDF1> bdf1_linsys(bdf1, geom.size_st());
  matvec_fd<BDF2> bdf2_linsys(bdf2, geom.size_st());

  newton_workspace nwsp(geom.size_st());
  GMRES_workspace  gwsp(settings.gmres, geom.size_st());

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

  write_hist(geom, hist, (nt + 1), dt, "hist.ks");

  delete[] state;
  delete[] hist;

  return 0;
}

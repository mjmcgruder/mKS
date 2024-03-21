#include "initialization.cpp"
#include "io.cpp"
#include "time.cpp"
#include "time_schemes.cpp"


int main(int argc, char** argv)
{
  newton_settings settings;
  settings.gmres.ninner = 100;
  settings.gmres.nouter = 5;
  settings.gmres.rtol   = 1e-3;
  settings.nonlin_tol   = 1e-8;
  settings.max_iter     = 10;

  double nc  = 1.;
  double nu2 = 1.;
  double nu4 = 0.1111;
  double eta = 16.;
  double sig = 4096.;
  double tau = 16.;

  // read_input("ifile");
  // exit(0);

  if (argc != 5)
  {
    printf("Usage: mKS <ne> <p> <nt> <dt>\n");
    exit(0);
  }

  const usize  ne = strtoull(argv[1], nullptr, 10);
  const usize  p  = strtoull(argv[2], nullptr, 10);
  const usize  nt = strtoull(argv[3], nullptr, 10);
  const double dt = strtod(argv[4], nullptr);

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

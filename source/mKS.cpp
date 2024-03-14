#include "initialization.cpp"
#include "io.cpp"
#include "time.cpp"
#include "time_schemes.cpp"


int main(int argc, char** argv)
{
  if (argc != 5)
  {
    printf("Usage: mKS <ne> <p> <nt> <dt>\n");
    exit(0);
  }

  const usize  ne = strtoull(argv[1], nullptr, 10);
  const usize  p  = strtoull(argv[2], nullptr, 10);
  const usize  nt = strtoull(argv[3], nullptr, 10);
  const double dt = strtod(argv[4], nullptr);

  comp_geom geom(1., ne, p);

  double* state = new double[geom.size_st()];
  double* hist  = new double[(nt + 1) * geom.size_st()];

  // init_constant(geom, 1.2, state);
  init_gaussian(geom, 1., 0.5, 0.1, state);

  // rk4(geom, nt, dt, state, hist);

  BDF1 bdf1(geom, dt);

  newton_settings settings{50, 10, 1e-8, 1e-8, 100};
  newton_workspace wsp(geom.size_st());

  bdf1.update_step(state);

  for (usize it = 0; it < nt; ++it)
  {
    printf("\n\n");
    printf("TIME STEP %zu\n", it);
    printf("\n\n");

    memcpy(hist + it * geom.size_st(), state, geom.size_st() * sizeof(double));

    newton_solve(bdf1, state, geom.size_st(), settings, wsp);

    bdf1.update_step(state);
  }

  memcpy(hist + nt * geom.size_st(), state, geom.size_st() * sizeof(double));

  write_hist(geom, hist, (nt + 1), "hist.ks");

  delete[] state;
  delete[] hist;

  return 0;
}

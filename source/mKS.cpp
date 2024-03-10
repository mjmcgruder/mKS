#include "time.cpp"
#include "initialization.cpp"
#include "io.cpp"


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

  rk4(geom, nt, dt, state, hist);

  write_hist(geom, hist, (nt + 1), "hist.ks");

  delete[] state;
  delete[] hist;

  return 0;
}

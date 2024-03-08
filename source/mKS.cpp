#include "time.cpp"
#include "initialization.cpp"
#include "io.cpp"


int main()
{
  const usize  ne = 16;
  const usize  p  = 2;
  const usize  nt = 20000;
  const double dt = 0.0001;

  comp_geom geom(1., ne, p);

  double* state = new double[geom.size_st()];
  double* hist  = new double[nt * geom.size_st()];

  // init_constant(geom, 1.2, state);
  init_gaussian(geom, 1., 0.5, 0.1, state);

  rk4(geom, nt, dt, state, hist);

  write_hist(geom, hist, nt, "hist.ks");

  delete[] state;
  delete[] hist;

  return 0;
}

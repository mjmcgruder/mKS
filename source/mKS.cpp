#include "time.cpp"
#include "initialization.cpp"
#include "io.cpp"


int main()
{
  const usize  ne = 10;
  const usize  p  = 2;
  const usize  nt = 10000;
  const double dt = 0.01;

  comp_geom geom(1., ne, p);

  double* state = new double[geom.size_st()];

  init_gaussian(geom, 1., 0.5, 0.1, state);

  rk4(geom, nt, dt, state);

  write_state(geom, state, "state.ks");

  delete[] state;

  return 0;
}

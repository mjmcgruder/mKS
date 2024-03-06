#include "time.cpp"
#include "initialization.cpp"
#include "io.cpp"


int main()
{
  const usize  ne = 100;
  const usize  p  = 7;
  const usize  nt = 1000000;
  const double dt = 0.000001;

  comp_geom geom(1., ne, p);

  double* state = new double[geom.size_st()];

  // init_constant(geom, 1.2, state);
  init_gaussian(geom, 1., 0.5, 0.1, state);

  rk4(geom, nt, dt, state);

  write_state(geom, state, "state.ks");

  delete[] state;

  return 0;
}

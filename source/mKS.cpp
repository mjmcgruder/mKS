#include "time.cpp"
#include "initialization.cpp"


int main()
{
  const usize  ne = 10;
  const usize  p  = 0;
  const usize  nt = 10;
  const double dt = 0.1;

  comp_geom geom(1., ne, p);

  double* state = new double[geom.size_st()];

  // init_constant(geom, 3.14, state);
  init_gaussian(geom, 1., 0.5, 0.2, state);

  euler(geom, nt, dt, state);

  printf("\n%f\n", state[0]);

  delete[] state;

  return 0;
}

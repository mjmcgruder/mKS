#include <cstdlib>

#include "test_quadrature.cpp"
#include "test_math.cpp"
#include "test_basis.cpp"
#include "test_geometry.cpp"
#include "test_linear_solver.cpp"
#include "test_newton.cpp"
#include "test_fourth.cpp"


#define RUN_TEST(func, name)         \
  {                                  \
    printf("TESTING %s...\n", name); \
    if (func != 0)                   \
    {                                \
      printf("  FAILED!\n");         \
    }                                \
    else                             \
    {                                \
      printf("  SUCCESS!\n");        \
    }                                \
  }


int main(int argc, char** argv)
{
  // RUN_TEST(test_quadrature(), "quadrature");
  // RUN_TEST(test_inverse(), "matrix inversion");
  // RUN_TEST(test_matmul(), "matrix multiplication");
  // RUN_TEST(test_basis(), "basis sum 1");
  // RUN_TEST(test_basis_d1(), "basis derivatives sum 0");
  // RUN_TEST(test_mass_matrix(), "mass matrix sum length");
  RUN_TEST(test_GMRES(), "GMRES (laplace)");
  RUN_TEST(test_newton(), "newton bowl");
  // RUN_TEST(test_convergence(), "fourth order term convergence")
}

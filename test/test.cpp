#include <cstdlib>

#include "test_quadrature.cpp"


#define RUN_TEST(func) \
  {                    \
    if (func != 0)     \
    {                  \
      exit(1);         \
    }                  \
  }


int main(int argc, char** argv)
{
  RUN_TEST(test_quadrature());
}

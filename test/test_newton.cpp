#pragma once


#include "time.cpp"


struct bowl
{
  // external facing vars

  double* Fo;

  // state info

  // --

  void uresid(double* u, double* F);
  void matvec(double* v, double* mv);

  void update_guess(double* x);

  // --

  bowl();
  ~bowl();

  bowl(const bowl& oth)            = delete;
  bowl(bowl&& oth)                 = delete;
  bowl& operator=(const bowl& oth) = delete;
  bowl& operator=(bowl&& oth)      = delete;
};


int test_newton()
{



  return 0;
}

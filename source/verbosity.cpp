#pragma once


#include <cstdio>


// 0 (or less): print nothing
// 1          : print almost nothing
// 2          : print a little
// 3 (or more): print a lot
static int VERBOSITY = 3;


#define vrbprintf(vrb, ...) \
  {                         \
    if ((vrb) <= VERBOSITY) \
    {                       \
      printf(__VA_ARGS__);  \
    }                       \
  }

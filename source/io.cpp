#pragma once


#include <cstdio>

#include "geometry.cpp"


void write_state(comp_geom& geom, double* state, const char* filename)
{
  FILE* f = fopen(filename, "w");

  fwrite(&geom.l,  sizeof(double), 1, f);
  fwrite(&geom.ne, sizeof(usize),  1, f);
  fwrite(&geom.p,  sizeof(usize),  1, f);

  fwrite(state, sizeof(double), geom.size_st(), f);

  fclose(f);
}


void write_hist(comp_geom& geom, double* hist, usize nt, double dt,
                const char* filename)
{
  FILE* f = fopen(filename, "w");

  fwrite(&geom.l,  sizeof(double), 1, f);
  fwrite(&geom.ne, sizeof(usize),  1, f);
  fwrite(&geom.p,  sizeof(usize),  1, f);
  fwrite(&nt,      sizeof(usize),  1, f);
  fwrite(&dt,      sizeof(double), 1, f);

  fwrite(hist, sizeof(double), nt * geom.size_st(), f);

  fclose(f);
}

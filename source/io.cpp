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


struct options
{
  // equation coefficients
  double a;
  double b;
  double c;
  // stabilization coefficients
  double eta;
  double sig;
  double tau;
  // temporal evolution
  usize  nt;
  double dt;
  // spatial resolution
  usize  ne;
  usize  p;
  double l;
  // nonlinear solver
  double nonlin_tol;
  double max_nonlin_iters;
  // linear solver
  usize  gmres_ninner;
  usize  gmres_nouter;
  double gmres_rtol;
};


#define LINE_BUF_SIZE 256


char peek(FILE* f)
{
  int ch = fgetc(f);
      ch = ungetc(ch, f);
  return (char)ch;
}


void seek_next_line(FILE* f)
{
  char c;
  do {
    c = fgetc(f);
  } while(c != '\n' && c != EOF);
}


bool read_keyval(char buf[LINE_BUF_SIZE],
                 char key[LINE_BUF_SIZE], char val[LINE_BUF_SIZE])
{
  int res = sscanf(buf, " %s : %s \n", key, val);

  bool read_ok = true;
  if (res != 2) { read_ok = false; }

  return read_ok;
}


// reads line up to comment (might return empty)
bool read_line(FILE* f, char buf[LINE_BUF_SIZE])
{
  char  c;
  usize bi = 0;
  for (; bi < LINE_BUF_SIZE - 1; ++bi)  // - 1 for null terminator
  {
    c = (char)fgetc(f);
    if (c == '#')       { seek_next_line(f); break; }
    else if (c == '\n') { break; }
    buf[bi] = c;
  }

  buf[bi] = '\0';
  ++bi;  // set final size

  bool read_ok = true;
  if (bi == LINE_BUF_SIZE)
  {
    char nextc = peek(f);
    if (nextc == '\n')     { c = (char)fgetc(f); }
    else if (nextc != EOF) { read_ok = false; }
  }

  return read_ok;
}


void read_input(const char* filename)
{
  FILE* f = fopen(filename, "r");
  if (f == nullptr)
  {
    printf("Failed to open input file %s\n", filename);
    exit(1);
  }

  usize line = 0;
  while((int)peek(f) != EOF)
  {
    ++line;

    char buf[LINE_BUF_SIZE];

    bool line_read_ok = read_line(f, buf);
    if (!line_read_ok)
    {
      printf("line %zu is over the %d character limit!\n",
             line, LINE_BUF_SIZE - 1);
      exit(1);
    }

    printf("%s\n", buf);

    if (strlen(buf) == 0) { continue; }

    char key[LINE_BUF_SIZE];
    char val[LINE_BUF_SIZE];
    bool keyval_read_ok = read_keyval(buf, key, val);
    if (!keyval_read_ok)
    {
      printf("Failed to read key value pair at line %zu\n", line);
      exit(1);
    }

    printf("%s\n", key);
    printf("%s\n", val);
  }

  fclose(f);
}

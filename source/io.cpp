#pragma once


#include <cstdio>
#include <map>


#include "geometry.cpp"


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


enum input_type
{
  itype_usize,
  itype_double
};


struct option
{
  input_type  type;
  void*       addr;
  const char* desc;
};


option mkopt(usize& data, const char* desc)
{
  return {itype_usize, &data, desc};
}

option mkopt(double& data, const char* desc)
{
  return {itype_double, &data, desc};
}


struct option_comp
{
  bool operator()(const char* a, const char* b) const
  {
    if (strcmp(a, b) < 0) { return true;  }
    else                  { return false; }
  }
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


void read_input(const char* filename,
                std::map<const char*, option, option_comp>& options)
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

    if (strlen(buf) == 0) { continue; }

    char key[LINE_BUF_SIZE];
    char val[LINE_BUF_SIZE];
    bool keyval_read_ok = read_keyval(buf, key, val);
    if (!keyval_read_ok)
    {
      printf("Failed to read key value pair at line %zu\n", line);
      exit(1);
    }

    auto it = options.find(key);
    if (it == options.end())
    {
      printf("Option \"%s\" not recognized\n", key);
      exit(1);
    }

    option optinfo = it->second;

    switch (optinfo.type)
    {
      case itype_usize:
        *((usize*)optinfo.addr) = (usize)strtoull(val, nullptr, 10);
        break;
      case itype_double:
        *((double*)optinfo.addr) = strtod(val, nullptr);
        break;
    }
  }

  fclose(f);
}

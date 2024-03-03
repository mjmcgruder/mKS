#pragma once


#include "geometry.cpp"


void compute_residual(comp_geom& geom, double* state, double* resid)
{
  // constants (set these somewhere else later)

  const double a = 1.;

  // zero residual

  memset(resid, 0, geom.size_st() * sizeof(double));

  // element interior flux

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    for (usize qi = 0; qi < geom.qrule.n; ++qi)
    {
      double s = 0.;
      for (usize bi = 0; bi < geom.p + 1; ++bi)
      {
        s += geom.iqp_bf(qi, bi) * geom.idx_st(state, ei, bi);
      }

      double flux = a * s;

      for (usize ti = 0; ti < geom.p + 1; ++ti)
      {
        geom.idx_st(resid, ei, ti) += geom.iqp_d1(qi, ti) * flux;
      }
    }
  }

  // element boundary flux (preiodic BCs)
  // each "element" deals with the boundary to its right

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    double el = ei;                  // left element
    double er = (ei + 1) % geom.ne;  // right element

    double sl = 0.;  // "left"  state, high side of left element
    double sr = 0.;  // "right" state, low side of right element

    for (usize bi = 0; bi < geom.p + 1; ++bi)
    {
      sl += geom.edg_bf(elem_side::high, bi) * geom.idx_st(state, el, bi);
      sr += geom.edg_bf(elem_side::low,  bi) * geom.idx_st(state, er, bi);
    }

    double fl   = a * sl;
    double fr   = a * sr;
    double flux = 0.5 * ((fl + fr) - std::abs(a) * (sr - sl));

    for (usize ti = 0; ti < geom.p + 1; ++ti)
    {
      geom.idx_st(resid, el, ti) += geom.edg_bf(elem_side::high, ti) * flux;
      geom.idx_st(resid, er, ti) -= geom.edg_bf(elem_side::low,  ti) * flux;
    }
  }
}

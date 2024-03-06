#pragma once


#include "geometry.cpp"


void compute_residual(comp_geom& geom, double* state, double* resid, double* f)
{
  // constants (set these somewhere else later)

  const double a = 1.;

  // zero residual

  memset(resid, 0, geom.size_st() * sizeof(double));

  // element interior flux

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double J   = geom.J;
    double iJ  = 1. / J;

    for (usize qi = 0; qi < geom.qrule.n; ++qi)
    {
      double s = 0.;
      for (usize bi = 0; bi < pp1; ++bi)
      {
        s += geom.iqp_bf(qi, bi) * geom.idx_st(state, ei, bi);
      }

      // // advective term
      // double advective_flux = a * s;

      // nonlinear term
      double nonlinear_flux = 0.5 * s * s;

      for (usize ti = 0; ti < pp1; ++ti)
      {
        // // advective term
        // geom.idx_st(resid, ei, ti) -=
        //   (geom.iqp_d1(qi, ti) * iJ) * advective_flux * J * geom.qrule.w[qi];

        // nonlinear term
        geom.idx_st(resid, ei, ti) -=
          (geom.iqp_d1(qi, ti) * iJ) * nonlinear_flux * J * geom.qrule.w[qi];
      }
    }
  }

  // element boundary flux (preiodic BCs)
  // each "element" deals with the boundary to its right

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double J   = geom.J;

    double el = ei;                  // left element
    double er = (ei + 1) % geom.ne;  // right element

    double sl = 0.;  // "left"  state, high side of left element
    double sr = 0.;  // "right" state, low side of right element

    for (usize bi = 0; bi < pp1; ++bi)
    {
      sl += geom.edg_bf(elem_side::high, bi) * geom.idx_st(state, el, bi);
      sr += geom.edg_bf(elem_side::low,  bi) * geom.idx_st(state, er, bi);
    }

    // // advective term
    // double ad_fl   = a * sl;
    // double ad_fr   = a * sr;
    // double ad_flux = 0.5 * ((ad_fl + ad_fr) - std::abs(a) * (sr - sl));

    // nonlinear term
    double nl_fl   = 0.5 * sl * sl;
    double nl_fr   = 0.5 * sr * sr;
    double nl_flux = 0.5 * ((nl_fl + nl_fr) - std::abs(a) * (sr - sl));

    for (usize ti = 0; ti < pp1; ++ti)
    {
      // // advective term
      // geom.idx_st(resid, el, ti) += geom.edg_bf(elem_side::high, ti) * ad_flux;
      // geom.idx_st(resid, er, ti) -= geom.edg_bf(elem_side::low,  ti) * ad_flux;

      // nonlinear term
      geom.idx_st(resid, el, ti) += geom.edg_bf(elem_side::high, ti) * nl_flux;
      geom.idx_st(resid, er, ti) -= geom.edg_bf(elem_side::low,  ti) * nl_flux;
    }
  }

  // include mass matrix

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double iJ  = 1. / geom.J;
    matmul(pp1, 1, pp1, -iJ, geom.Minv, resid + ei * pp1, 0., f + ei * pp1);
  }
}

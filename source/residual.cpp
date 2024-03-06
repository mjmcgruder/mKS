#pragma once


#include "geometry.cpp"


void compute_residual(comp_geom& geom, double* state, double* resid, double* f)
{
  // constants (set these somewhere else later)

  const double a   = 1.;
  const double nu  = 0.1;
  const double eta = 32.;

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
      double s   = 0.;
      double s_x = 0.;
      for (usize bi = 0; bi < pp1; ++bi)
      {
        s   += geom.iqp_bf(qi, bi) * geom.idx_st(state, ei, bi);

        s_x += (geom.iqp_d1(qi, bi) * iJ) * geom.idx_st(state, ei, bi);
      }

      // // for advective term
      // double advective_flux = a * s;

      // // for nonlinear term
      // double nonlinear_flux = 0.5 * s * s;

      // for viscous term
      double viscous_flux = nu * s_x;

      for (usize ti = 0; ti < pp1; ++ti)
      {
        // // advective term
        // geom.idx_st(resid, ei, ti) -=
        //   (geom.iqp_d1(qi, ti) * iJ) * advective_flux * J * geom.qrule.w[qi];

        // // nonlinear term
        // geom.idx_st(resid, ei, ti) -=
        //   (geom.iqp_d1(qi, ti) * iJ) * nonlinear_flux * J * geom.qrule.w[qi];

        // viscous term
        geom.idx_st(resid, ei, ti) +=
          (geom.iqp_d1(qi, ti) * iJ) * viscous_flux * J * geom.qrule.w[qi];

      }
    }
  }

  // element boundary flux (preiodic BCs)
  // each "element" deals with the boundary to its right

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double J   = geom.J;
    double iJ  = 1. / J;

    double el = ei;                  // left element
    double er = (ei + 1) % geom.ne;  // right element (with wrap around)

    double sl   = 0.;  // "left"  state, high side of left element
    double sr   = 0.;  // "right" state, low side of right element
    double sl_x = 0.;
    double sr_x = 0.;

    for (usize bi = 0; bi < pp1; ++bi)
    {
      sl   += geom.edg_bf(elem_side::high, bi) * geom.idx_st(state, el, bi);
      sr   += geom.edg_bf(elem_side::low,  bi) * geom.idx_st(state, er, bi);

      sl_x += (geom.edg_d1(elem_side::high, bi) * iJ) * geom.idx_st(state, el, bi);
      sr_x += (geom.edg_d1(elem_side::low,  bi) * iJ) * geom.idx_st(state, er, bi);
    }

    // // for advective term
    // double ad_fl   = a * sl;
    // double ad_fr   = a * sr;
    // double ad_flux = 0.5 * ((ad_fl + ad_fr) - std::abs(a) * (sr - sl));

    // // for nonlinear term
    // double nl_fl   = 0.5 * sl * sl;
    // double nl_fr   = 0.5 * sr * sr;
    // double nl_flux = 0.5 * ((nl_fl + nl_fr) - std::abs(a) * (sr - sl));

    // for dual consistency term
    double s_av  = 0.5 * (sl + sr);
    double dc_fl = nu * (sl - s_av) * +1.;
    double dc_fr = nu * (sr - s_av) * -1.;

    // for viscous flux term
    double dl = (1. / geom.elem_len(el)) * nu * (sl - sr) * +1.;
    double dr = (1. / geom.elem_len(er)) * nu * (sr - sl) * -1.;
    double dm = 0.5 * (dl + dr);
    double ql = nu * sl_x;
    double qr = nu * sr_x;
    double qm = 0.5 * (ql + qr);
    double qh = qm - eta * dm;

    for (usize ti = 0; ti < pp1; ++ti)
    {
      // // advective term
      // geom.idx_st(resid, el, ti) += geom.edg_bf(elem_side::high, ti) * ad_flux;
      // geom.idx_st(resid, er, ti) -= geom.edg_bf(elem_side::low,  ti) * ad_flux;

      // // nonlinear term
      // geom.idx_st(resid, el, ti) += geom.edg_bf(elem_side::high, ti) * nl_flux;
      // geom.idx_st(resid, er, ti) -= geom.edg_bf(elem_side::low,  ti) * nl_flux;

      // dual consistency term
      geom.idx_st(resid, el, ti) -= geom.edg_d1(elem_side::high, ti) * iJ * dc_fl;
      geom.idx_st(resid, er, ti) -= geom.edg_d1(elem_side::low,  ti) * iJ * dc_fr;

      // viscous flux term
      geom.idx_st(resid, el, ti) -= geom.edg_bf(elem_side::high, ti) * qh * +1.;
      geom.idx_st(resid, er, ti) -= geom.edg_bf(elem_side::low,  ti) * qh * -1.;
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

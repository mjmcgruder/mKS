#pragma once


#include "geometry.cpp"


void compute_residual(comp_geom& geom, double* state, double* resid, double* f)
{
  // zero residual

  memset(resid, 0, geom.size_st() * sizeof(double));

  // element interior flux

  for (usize ei = 0; ei < geom.ne; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double J   = geom.J;
    double iJ  = 1. / J;
    double iJ2 = iJ * iJ;

    for (usize qi = 0; qi < geom.qrule.n; ++qi)
    {
      double s    = 0.;
      double s_x  = 0.;
      double s_xx = 0.;
      for (usize bi = 0; bi < pp1; ++bi)
      {
        double U = geom.idx_st(state, ei, bi);
        s       += (geom.iqp_bf(qi, bi)      ) * U;
        s_x     += (geom.iqp_d1(qi, bi) * iJ ) * U;
        s_xx    += (geom.iqp_d2(qi, bi) * iJ2) * U;
      }

      // for nonlinear term
      double nonlinear_flux = geom.nc * 0.5 * s * s;
      // for viscous term
      double viscous_flux = geom.nu2 * s_x;
      // for fourth order term
      double fourth_flux = geom.nu4 * s_xx;

      for (usize ti = 0; ti < pp1; ++ti)
      {
        double t_x  = geom.iqp_d1(qi, ti) * iJ;
        double t_xx = geom.iqp_d2(qi, ti) * iJ2;

        // nonlinear term
        geom.idx_st(resid, ei, ti) -= t_x * nonlinear_flux * J * geom.qrule.w[qi];
        // second order viscous term
        geom.idx_st(resid, ei, ti) -= t_x * viscous_flux * J * geom.qrule.w[qi];
        // fourth order viscous term
        geom.idx_st(resid, ei, ti) += t_xx * fourth_flux * J * geom.qrule.w[qi];

      }
    }
  }

  // source term

  if (geom.source != nullptr)
  {
    for (usize ei = 0; ei < geom.ne; ++ei)
    {
      usize  pp1 = geom.p + 1;
      double J   = geom.J;

      for (usize qi = 0; qi < geom.qrule.n; ++qi)
      {
        double x = geom.ref2glo(ei, geom.qrule.x[qi]);
        double s = geom.source(x);

        for (usize ti = 0; ti < pp1; ++ti)
        {
          geom.idx_st(resid, ei, ti) +=
            geom.iqp_bf(qi, ti) * s * J * geom.qrule.w[qi];
        }
      }
    }
  }

  // element boundary flux

  usize start_edge = 1;
  if (geom.btype == dirichlet) { start_edge = 0; }

  for (usize ei = start_edge; ei < geom.ne + 1; ++ei)
  {
    usize  pp1 = geom.p + 1;
    double J   = geom.J;
    double iJ  = 1. / J;
    double iJ2 = iJ * iJ;
    double iJ3 = iJ * iJ * iJ;

    int el = ei - 1;        // left element
    int er = ei % geom.ne;  // right element (with wrap around)
    if (ei == 0 && (geom.btype == dirichlet))       { el = -1; }
    if (ei == geom.ne && (geom.btype == dirichlet)) { er = -1; }

    double sl = 0.;  // "left"  state, high side of left element
    double sr = 0.;  // "right" state, low side of right element
    double sl_x   = 0.; double sr_x   = 0.;
    double sl_xx  = 0.; double sr_xx  = 0.;
    double sl_xxx = 0.; double sr_xxx = 0.;

    for (usize bi = 0; bi < pp1; ++bi)
    {
      if (el >= 0)
      {
        double Ul = geom.idx_st(state, el, bi);
        sl       += (geom.edg_bf(elem_side::high, bi)      ) * Ul;
        sl_x     += (geom.edg_d1(elem_side::high, bi) * iJ ) * Ul;
        sl_xx    += (geom.edg_d2(elem_side::high, bi) * iJ2) * Ul;
        sl_xxx   += (geom.edg_d3(elem_side::high, bi) * iJ3) * Ul;
      }

      if (er >= 0)
      {
        double Ur = geom.idx_st(state, er, bi);
        sr       += (geom.edg_bf(elem_side::low,  bi)      ) * Ur;
        sr_x     += (geom.edg_d1(elem_side::low,  bi) * iJ ) * Ur;
        sr_xx    += (geom.edg_d2(elem_side::low,  bi) * iJ2) * Ur;
        sr_xxx   += (geom.edg_d3(elem_side::low,  bi) * iJ3) * Ur;
      }
    }

    if (el < 0)
    {
      sl     = geom.value_l;
      sl_x   = geom.deriv_l;
      sl_xx  = sr_xx;
      sl_xxx = sr_xxx;
    }

    if (er < 0)
    {
      sr     = geom.value_r;
      sr_x   = geom.deriv_r;
      sr_xx  = sl_xx;
      sr_xxx = sl_xxx;
    }

    double ihl = 1. / geom.elem_len(el);
    double ihr = 1. / geom.elem_len(er);

    // for nonlinear term
    double nl_fl   = 0.5 * sl * sl;
    double nl_fr   = 0.5 * sr * sr;
    double ah;
    if (sl == sr) { ah = sl; }
    else          { ah = (nl_fr - nl_fl) / (sr - sl); }
    double nl_flux = geom.nc * 0.5 * ((nl_fl + nl_fr) - std::abs(ah) * (sr - sl));

    // for second order viscous term
    double s_av  = 0.5 * (sl + sr);
    double dc_fl = geom.nu2 * (sl - s_av) * +1.;
    double dc_fr = geom.nu2 * (sr - s_av) * -1.;
    double dl    = ihl * geom.nu2 * (sl - sr) * +1.;
    double dr    = ihr * geom.nu2 * (sr - sl) * -1.;
    double dm    = 0.5 * (dl + dr);
    double ql    = geom.nu2 * sl_x;
    double qr    = geom.nu2 * sr_x;
    double qm    = 0.5 * (ql + qr);
    double qh    = qm - geom.eta * dm;

    // for fourth order viscous term
    double avg_s_xx  = 0.5 * (sl_xx + sr_xx);
    double avg_s_xxx = 0.5 * (sl_xxx + sr_xxx);
    double diff_sl    = sl - sr;
    double diff_sr    = sr - sl;
    double diff_sl_x  = sl_x - sr_x;
    double diff_sr_x  = sr_x - sl_x;
    double sigl       = geom.sig * ihl * ihl * ihl;
    double taul       = geom.tau * ihl;
    double sigr       = geom.sig * ihr * ihr * ihr;
    double taur       = geom.tau * ihr;

    for (usize ti = 0; ti < pp1; ++ti)
    {
      if (el >= 0)
      {
        double tl     = geom.edg_bf(elem_side::high, ti);
        double tl_x   = geom.edg_d1(elem_side::high, ti) * iJ;
        double tl_xx  = geom.edg_d2(elem_side::high, ti) * iJ2;
        double tl_xxx = geom.edg_d3(elem_side::high, ti) * iJ3;
        // nonlinear term
        geom.idx_st(resid, el, ti) += tl * nl_flux;
        // second order viscous term
        geom.idx_st(resid, el, ti) += tl_x * dc_fl;
        geom.idx_st(resid, el, ti) += tl * qh * +1.;
        // fourth order viscous term
        geom.idx_st(resid, el, ti) +=
          geom.nu4 *
            (avg_s_xxx * tl * +1 +
             0.5 * tl_xxx * diff_sl * +1 -
             avg_s_xx * tl_x * +1 -
             0.5 * tl_xx * diff_sl_x * +1 +
             sigl * diff_sl * tl * +1 +
             taul * diff_sl_x * tl_x * +1);
      }

      if (er >= 0)
      {
        double tr     = geom.edg_bf(elem_side::low,  ti);
        double tr_x   = geom.edg_d1(elem_side::low,  ti) * iJ;
        double tr_xx  = geom.edg_d2(elem_side::low,  ti) * iJ2;
        double tr_xxx = geom.edg_d3(elem_side::low,  ti) * iJ3;
        // nonlinear term
        geom.idx_st(resid, er, ti) -= tr * nl_flux;
        // second order viscous term
        geom.idx_st(resid, er, ti) += tr_x * dc_fr;
        geom.idx_st(resid, er, ti) += tr * qh * -1.;
        // fourth order viscous term
        geom.idx_st(resid, er, ti) +=
          geom.nu4 *
            (avg_s_xxx * tr * -1 +
             0.5 * tr_xxx * diff_sr * -1 -
             avg_s_xx * tr_x * -1 -
             0.5 * tr_xx * diff_sr_x * -1 +
             sigr * diff_sr * tr * +1 +
             taur * diff_sr_x * tr_x * +1);
      }
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

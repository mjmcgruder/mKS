#pragma once


#include "basis.cpp"
#include "linalg.cpp"
#include "quadrature.cpp"


enum struct elem_side
{
  low  = 0,
  high = 1,
};


// Structure containing discretization information with the resulting
// pre-computed quantities. For now, all elements are the same size and the same
// order to keep things simple.
struct comp_geom
{
  // KS equation parameters (no equation set distinction in this code) ---

  double nc;   // nonlinear coefficient
  double nu2;  // second order viscosity coefficient
  double nu4;  // fourth order viscosity coefficient

  double eta;  // second order IP stabilization constant
  double sig;  // fourth order IP stabilization constants
  double tau;  // *

  // geometry ---

  double l;   // domain length
  usize  ne;  // number of elements
  usize  p;   // element order (all elements same order for now)

  // Jacobian / determinant
  // Only one for now, all elements are same size and geometry within elements
  // is linear. Change this to a per-element Jacobian later.
  double J;

  double elem_len(usize e);
  double ref2glo(usize e, double xi);

  // pre-computes ---

  // Quadrature rule sufficient for integrating polynomials of order 2 * p + 1
  // exactly. This should turn into a set of different quadrature rules when
  // variable order is supported. Or maybe just a single highest-order rule.
  quadrature_rule qrule;

  // mass matrix
  // Only a single mass matrix for now since all elements are the same order and
  // the Jacobian determinant is constant wihtin elements. Change this to one
  // mass matrix per order later.
  double* Minv;

  // pre-computed element interior basis functions and derivatives [qp [bf]]
  double* iqp_bf_arr;
  double* iqp_d1_arr;
  double* iqp_d2_arr;

  // pre-computed element edge basis functions and derivatives [l/h [bf]]
  double* edg_bf_arr;
  double* edg_d1_arr;
  double* edg_d2_arr;
  double* edg_d3_arr;

  // array indexing ---

  double& iqp_bf(usize qp, usize bf);
  double& iqp_d1(usize qp, usize bf);
  double& iqp_d2(usize qp, usize bf);
  double& edg_bf(elem_side side, usize bf);
  double& edg_d1(elem_side side, usize bf);
  double& edg_d2(elem_side side, usize bf);
  double& edg_d3(elem_side side, usize bf);

  double& idx_st(double* state, usize e, usize bf);

  // array sizing ---

  usize size_st();

  // constructors etc. ---

  comp_geom(double l_, usize ne_, usize p_, double nc_, double nu2_,
            double nu4_, double eta_, double sig_, double tau_);
  ~comp_geom();

  comp_geom(const comp_geom& oth)            = delete;
  comp_geom(comp_geom&& oth)                 = delete;
  comp_geom& operator=(const comp_geom& oth) = delete;
  comp_geom& operator=(comp_geom&& oth)      = delete;
};


double comp_geom::elem_len(usize e)
{
  return l / double(ne);
}

double comp_geom::ref2glo(usize e, double xi)
{
  double elen = elem_len(e);
  return (double(e) * elen) + ((xi + 1.) * J);
}


double& comp_geom::iqp_bf(usize qp, usize bf)
{
  return iqp_bf_arr[(p + 1) * qp + bf];
}

double& comp_geom::iqp_d1(usize qp, usize bf)
{
  return iqp_d1_arr[(p + 1) * qp + bf];
}

double& comp_geom::iqp_d2(usize qp, usize bf)
{
  return iqp_d2_arr[(p + 1) * qp + bf];
}

double& comp_geom::edg_bf(elem_side side, usize bf)
{
  return edg_bf_arr[(p + 1) * usize(side) + bf];
}

double& comp_geom::edg_d1(elem_side side, usize bf)
{
  return edg_d1_arr[(p + 1) * usize(side) + bf];
}

double& comp_geom::edg_d2(elem_side side, usize bf)
{
  return edg_d2_arr[(p + 1) * usize(side) + bf];
}

double& comp_geom::edg_d3(elem_side side, usize bf)
{
  return edg_d3_arr[(p + 1) * usize(side) + bf];
}

double& comp_geom::idx_st(double* state, usize e, usize bf)
{
  return state[(p + 1) * e + bf];
}


usize comp_geom::size_st()
{
  return ne * (p + 1);
}


comp_geom::comp_geom(double l_, usize ne_, usize p_, double nc_, double nu2_,
                     double nu4_, double eta_, double sig_, double tau_)
{
  nc  = nc_;
  nu2 = nu2_;
  nu4 = nu4_;
  eta = eta_;
  sig = sig_;
  tau = tau_;

  l  = l_;
  ne = ne_;
  p  = p_;

  qrule = get_quadrature_rule(p);

  J = 0.5 * elem_len(0);

  // interior basis functions and gradients

  iqp_bf_arr = new double[qrule.n * (p + 1)];
  iqp_d1_arr = new double[qrule.n * (p + 1)];
  iqp_d2_arr = new double[qrule.n * (p + 1)];

  for (usize iq = 0; iq < qrule.n; ++iq)
  {
    lagrange_chebyshev_bf(p, qrule.x[iq], &iqp_bf(iq, 0));
    lagrange_chebyshev_d1(p, qrule.x[iq], &iqp_d1(iq, 0));
    lagrange_chebyshev_d2(p, qrule.x[iq], &iqp_d2(iq, 0));
  }

  // edge basis functions and gradients

  edg_bf_arr = new double[2 * (p + 1)];
  edg_d1_arr = new double[2 * (p + 1)];
  edg_d2_arr = new double[2 * (p + 1)];
  edg_d3_arr = new double[2 * (p + 1)];

  lagrange_chebyshev_bf(p, -1., &edg_bf(elem_side::low, 0));
  lagrange_chebyshev_bf(p, +1., &edg_bf(elem_side::high, 0));

  lagrange_chebyshev_d1(p, -1., &edg_d1(elem_side::low, 0));
  lagrange_chebyshev_d1(p, +1., &edg_d1(elem_side::high, 0));

  lagrange_chebyshev_d2(p, -1., &edg_d2(elem_side::low, 0));
  lagrange_chebyshev_d2(p, +1., &edg_d2(elem_side::high, 0));

  lagrange_chebyshev_d3(p, -1., &edg_d3(elem_side::low, 0));
  lagrange_chebyshev_d3(p, +1., &edg_d3(elem_side::high, 0));

  // make and invert mass matrix

  Minv      = new double[(p + 1) * (p + 1)];
  double* M = new double[(p + 1) * (p + 1)];

  for (usize tf = 0; tf < p + 1; ++tf)
  {
    for (usize bf = tf; bf < p + 1; ++bf)
    {
      double acc = 0.;
      for (usize iq = 0; iq < qrule.n; ++iq)
      {
        acc += iqp_bf(iq, tf) * iqp_bf(iq, bf) * qrule.w[iq];
      }
      M[(p + 1) * tf + bf] = acc;
      M[(p + 1) * bf + tf] = acc;
    }
  }

  invert_matrix(M, p + 1, Minv);

  delete[] M;
}


comp_geom::~comp_geom()
{
  delete[] Minv;
  delete[] iqp_bf_arr;
  delete[] iqp_d1_arr;
  delete[] iqp_d2_arr;
  delete[] edg_bf_arr;
  delete[] edg_d1_arr;
  delete[] edg_d2_arr;
  delete[] edg_d3_arr;
}

#pragma once


#include "types.cpp"


// // TODO: in future, template this and provide time scheme specific
// // specializations, maybe a value template with time scheme enum or something
// struct BDF1
// {
//   // simulation information
//   comp_geom& geom;
//   double     dt;
//
//   // known time steps (owned)
//   double* un;
//
//   // linearization state (current nonlinear solve guess) (not owned)
//   double* uk;
//
//   // residual work space (owned)
//   double* resid;
//   double* f;        // steady residual
//   double* uk_pert;  // perturbed linearization state, for finite difference
//
//   double* Fo;       // original (linearization state) unsteady residual
//   double* Fp;       // perturbed state unsteady residual
//
//   // computational routines
//
//   void uresid(double* u, double* F);   // (unsteady) resiudal
//   void matvec(double* v, double* mv);  // jacobian vector product
//
//   // internal state update routines
//
//   void update_linearization_state(double* st);
//   void update_time_step(double* st);
//
//   // --
//
//   BDF1(comp_geom& geom, double dt);
//   ~BDF1();
//
//   BDF1(const BDF1& oth)            = delete;
//   BDF1(BDF1&& oth)                 = delete;
//   BDF1& operator=(const BDF1& oth) = delete;
//   BDF1& operator=(BDF1&& oth)      = delete;
// };
//
//
// void BDF1::uresid(double* u, double* F)
// {
//   compute_residual(geom, u, resid, f);
//
//   for (usize i = 0; i < geom.size_st(); ++i)
//   {
//     F[i] = u[i] - un[i] - dt * f[i];
//   }
// }
//
//
// void BDF1::matvec(double* v, double* mv)
// {
//   double normv = norm(v, geom.size_st()) + DBL_EPSILON;
//   double eps   = 1e-4;
//
//   for (usize i = 0; i < geom.size_st(); ++i)
//   {
//     uk_pert[i] = uk[i] + eps * (v[i] / normv);
//   }
//
//   uresid(uk_pert, Fp);
//
//   for (usize i = 0; i < geom.size_st(); ++i)
//   {
//     mv[i] = (Fp[i] - Fo[i]) / eps;
//   }
// }
//
//
// void BDF1::update_linearization_state(double* new_uk)
// {
//   uk = new_uk;
//   uresid(uk, Fo);
// }
//
//
// void BDF1::update_time_step(double* u)
// {
//   memcpy(un, u, geom.size_st() * sizeof(*u));
// }
//
//
// BDF1::BDF1(comp_geom& geom_, double dt_) : geom(geom_), dt(dt_)
// {
//   un      = new double[geom.size_st()];
//   resid   = new double[geom.size_st()];
//   f       = new double[geom.size_st()];
//   uk_pert = new double[geom.size_st()];
//   Fo      = new double[geom.size_st()];
//   Fp      = new double[geom.size_st()];
// }
//
//
// BDF1::~BDF1()
// {
//   delete[] un;
//   delete[] resid;
//   delete[] f;
//   delete[] uk_pert;
//   delete[] Fo;
//   delete[] Fp;
// }

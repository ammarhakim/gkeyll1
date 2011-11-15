/**
 * @file	LcTenMomentEquation.cpp
 *
 * @brief	TenMoment equations for gas-dynamics.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTenMomentEquation.h>

// std includes
#include <cmath>

namespace Lucee
{
// set ids for creators
  const char *TenMomentEquation::id = "TenMoment";

// makes indexeing easier
  static const unsigned RHO = 0;
  static const unsigned U1 = 1;
  static const unsigned U2 = 2;
  static const unsigned U3 = 3;
  static const unsigned P11 = 4;
  static const unsigned P12 = 5;
  static const unsigned P13 = 6;
  static const unsigned P22 = 7;
  static const unsigned P23 = 8;
  static const unsigned P33 = 9;

  TenMomentEquation::TenMomentEquation()
    : Lucee::HyperEquation(10, 5)
  {
// equation system has 10 equations and 5 waves
  }
  
  void
  TenMomentEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
  }

  void
  TenMomentEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToLocal(&inQ[1], &outQ[1]); // momentum
// TODO PRESSURE TENSOR
  }

  void
  TenMomentEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToGlobal(&inQ[1], &outQ[1]); // momentum
// TODO PRESSURE TENSOR
  }

  void
  TenMomentEquation::flux(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
    Lucee::FieldPtr<double> v(10);
// compute primitive variables first
    primitive(q, v);
// density flux
    f[RHO] = v[RHO]*v[U1];
// momentum density flux
    f[U1] = v[RHO]*v[U1]*v[U1] + v[P11];
    f[U2] = v[RHO]*v[U1]*v[U2] + v[P12];
    f[U3] = v[RHO]*v[U1]*v[U3] + v[P13];
// total pressure flux
    f[P11] = v[RHO]*v[U1]*v[U1]*v[U1] + 3*v[U1]*v[P13];
    f[P12] = v[RHO]*v[U1]*v[U1]*v[U2] + 2*v[U1]*v[P12] + v[U2]*v[P11];
    f[P13] = v[RHO]*v[U1]*v[U1]*v[U3] + 2*v[U1]*v[P13] + v[U3]*v[P11];
    f[P22] = v[RHO]*v[U1]*v[U2]*v[U2] + v[U1]*v[P22] + 2*v[U2]*v[P12];
    f[P23] = v[RHO]*v[U1]*v[U2]*v[U3] + v[U1]*v[P23] + v[U2]*v[P13] + v[U3]*v[P12];
    f[P33] = v[RHO]*v[U1]*v[U3]*v[U3] + v[U1]*v[P33] + 2*v[U3]*v[P13];
  }

  void
  TenMomentEquation::speeds(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, double s[2])
  {
// TODO
  }

  void
  TenMomentEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v) const
  {
    double rho = q[RHO];
    v[RHO] = rho; // rho
// velocity components
    v[U1] = q[U1]/rho;
    v[U2] = q[U2]/rho;
    v[U3] = q[U3]/rho;
// pressure tensor
    v[P11] = q[P11] - rho*v[U1]*v[U1];
    v[P12] = q[P12] - rho*v[U1]*v[U2];
    v[P13] = q[P13] - rho*v[U1]*v[U3];
    v[P22] = q[P22] - rho*v[U2]*v[U2];
    v[P23] = q[P23] - rho*v[U2]*v[U3];
    v[P33] = q[P33] - rho*v[U3]*v[U3];
  }

  void
  TenMomentEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q) const
  {
    double rho = v[RHO];
    q[RHO] = rho; // rho
// velocity components
    q[U1] = rho*v[U1];
    q[U2] = rho*v[U2];
    q[U3] = rho*v[U3];
// pressure tensor
    q[P11] = v[P11] + rho*v[U1]*v[U1];
    q[P12] = v[P12] + rho*v[U1]*v[U2];
    q[P13] = v[P13] + rho*v[U1]*v[U3];
    q[P22] = v[P22] + rho*v[U2]*v[U2];
    q[P23] = v[P23] + rho*v[U2]*v[U3];
    q[P33] = v[P33] + rho*v[U3]*v[U3];
  }

  void
  TenMomentEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
// TODO
  }
}

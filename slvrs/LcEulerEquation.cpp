/**
 * @file	LcEulerEquation.cpp
 *
 * @brief	Euler equations for gas-dynamics.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEulerEquation.h>

// std includes
#include <cmath>

namespace Lucee
{
// set ids for creators
  const char *EulerEquation::id = "Euler";

  EulerEquation::EulerEquation()
    : Lucee::HyperEquation(5, 3)
  {
  }
  
  void
  EulerEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
    if (tbl.hasNumber("gasGamma"))
      gas_gamma = tbl.getNumber("gasGamma");
    else
      throw Lucee::Except("EulerEquation::readInput: Must specify gasGamma (gas constant)");
  }

  void
  EulerEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToLocal(&inQ[1], &outQ[1]); // momentum
    outQ[4] = inQ[4]; // energy density
  }

  void
  EulerEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToGlobal(&inQ[1], &outQ[1]); // momentum
    outQ[4] = inQ[4]; // energy density
  }

  void
  EulerEquation::flux(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
// compute pressure
    double pr = pressure(q);
// now compute flux
    f[0] = q[1]; // rho*u
    f[1] = q[1]*q[1]/q[0] + pr; // rho*u*u + pr
    f[2] = q[1]*q[2]/q[0]; // rho*u*v
    f[3] = q[1]*q[3]/q[0]; // rho*u*w
    f[4] = (q[4]+pr)*q[1]/q[0]; // (E+pr)*u
  }

  void
  EulerEquation::speeds(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, double s[2])
  {
// compute pressure
    double pr = pressure(q);
    double cs = std::sqrt(gas_gamma*pr/q[0]); // sound speed
    double u = q[1]/q[0]; // fluid velocity
    s[0] = u-cs;
    s[1] = u+cs;
  }

  void
  EulerEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v) const
  {
    v[0] = q[0]; // rho
    v[1] = q[1]/q[0]; // u
    v[2] = q[2]/q[0]; // v
    v[3] = q[3]/q[0]; // w
    v[4] = pressure(q); // pressure
  }

  void
  EulerEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q) const
  {
    q[0] = v[0]; // rho
    q[1] = v[0]*v[1]; // rho*u
    q[2] = v[0]*v[2]; // rho*v
    q[3] = v[0]*v[3]; // rho*w
    q[4] = v[4]/(gas_gamma-1) + 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]); // energy
  }

  void
  EulerEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double gas_gamma1 = gas_gamma-1;

// compute Roe averages for use in Riemann solver
    double rhsqrtl = std::sqrt(ql[0]); // left sqrt(rho)
    double rhsqrtr = std::sqrt(qr[0]); // right sqrt(rho)
    double pl = pressure(ql); // left pressure
    double pr = pressure(qr); // right pressure

    double rhsq2 = rhsqrtl + rhsqrtr;
// compute Roe averaged velocity components
    double u = (ql[1]/rhsqrtl  + qr[1]/rhsqrtr)/rhsq2;
    double v = (ql[2]/rhsqrtl  + qr[2]/rhsqrtr)/rhsq2;
    double w = (ql[3]/rhsqrtl  + qr[3]/rhsqrtr)/rhsq2;

    double q2 = u*u + v*v + w*w;
// Roe averaged enthalpy
    double enth = ((ql[4]+pl)/rhsqrtl + (qr[4]+pr)/rhsqrtr)/rhsq2;
// speed of sound
    double aa2 = gas_gamma1*(enth - 0.5*q2);
    double a = std::sqrt(aa2);
// other quantities
    double g1a2 = gas_gamma1/aa2;
    double euv = enth - q2;

// project onto left eigenvectors with Roe averaged values (see Tech Note 1007)
    double a3 = g1a2*(euv*jump[0] + u*jump[1] + v*jump[2] + w*jump[3] - jump[4]);
    double a1 = -v*jump[0] + jump[2];
    double a2 = -w*jump[0] + jump[3];
    double a4 = (jump[1] + (a-u)*jump[0] - a*a3)/(2*a);
    double a0 = jump[0] - a3 - a4;

// compute waves (see Tech Note 1007)

// wave 1: eigenvalue is u-c
    waves(0,0) = a0;
    waves(1,0) = a0*(u-a);
    waves(2,0) = a0*v;
    waves(3,0) = a0*w;
    waves(4,0) = a0*(enth-u*a);
    s[0] = u-a;

// wave 2: 3 eigenvectors for the repeated roots are combined into a
// single wave
    waves(0,1) = a3;
    waves(1,1) = a3*u;
    waves(2,1) = a1 + v*a3;
    waves(3,1) = a2 + w*a3;
    waves(4,1) = a1*v + a2*w + 0.5*a3*q2;
    s[1] = u;

// wave 3: eigenvalue is u+c
    waves(0,2) = a4;
    waves(1,2) = a4*(u+a);
    waves(2,2) = a4*v;
    waves(3,2) = a4*w;
    waves(4,2) = a4*(enth + u*a);
    s[2] = u+a;
  }

  double
  EulerEquation::pressure(const Lucee::ConstFieldPtr<double>& q) const
  {
    return (gas_gamma-1)*(q[4] - 0.5*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3])/q[0]);
  }
}

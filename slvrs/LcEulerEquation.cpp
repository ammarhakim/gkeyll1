/**
 * @file	LcEulerEquation.cpp
 *
 * @brief	Euler equations for gas-dynamics.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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
    if (tbl.hasNumber("gas_gamma"))
      gas_gamma = tbl.getNumber("gas_gamma");
    else
      throw Lucee::Except("EulerEquation::readInput: Must specify gas_gamma (gas constant)");
  }

  void
  EulerEquation::rotateToLocal(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& inQ, Lucee::FieldPtr<double>& outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToLocal(&inQ[1], &outQ[1]); // momentum
    outQ[4] = inQ[4]; // energy density
  }

  void
  EulerEquation::rotateToGlobal(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& inQ, Lucee::FieldPtr<double>& outQ)
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
  EulerEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v)
  {
    v[0] = q[0]; // rho
    v[1] = q[1]/q[0]; // u
    v[2] = q[2]/q[0]; // v
    v[3] = q[3]/q[0]; // w
    v[4] = pressure(q); // pressure
  }

  void
  EulerEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q)
  {
    q[0] = v[0]; // rho
    q[1] = v[0]*v[1]; // rho*u
    q[2] = v[0]*v[2]; // rho*v
    q[3] = v[0]*v[3]; // rho*w
    q[4] = v[4]/(gas_gamma-1) + 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]); // pressure
  }

  double
  EulerEquation::pressure(const Lucee::ConstFieldPtr<double>& q) const
  {
    return (gas_gamma-1)*(q[4] - 0.5*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3])/q[0]);
  }
}

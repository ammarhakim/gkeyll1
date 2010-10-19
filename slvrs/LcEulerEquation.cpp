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
    gas_gamma = tbl.getNumber("gas_gamma");
  }

  void
  EulerEquation::flux(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
// compute pressure
    double pr = (gas_gamma-1)*(q[4] - 0.5*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3])/q[0]);
// compute flux
    f[0] = q[1]; // rho*u
    f[1] = q[1]*q[1]/q[0] + pr; // rho*u*u + pr
    f[2] = q[1]*q[2]/q[0]; // rho*u*v
    f[3] = q[1]*q[3]/q[0]; // rho*u*w
    f[4] = (q[4]+pr)*q[1]/q[0]; // (E+pr)*u
  }
}

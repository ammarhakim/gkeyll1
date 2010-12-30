/**
 * @file	LcAdvectionEquation.cpp
 *
 * @brief	Advection equation
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
#include <LcAdvectionEquation.h>

// std includes
#include <algorithm>

namespace Lucee
{
  const char *AdvectionEquation::id = "Advection";

  AdvectionEquation::AdvectionEquation()
    : Lucee::HyperEquation(1, 1)
  {
  }

  void
  AdvectionEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
    std::vector<double> s;
    u[0] = u[1] = u[2] = 0.0; // by default set speeds to 0.0
    if (tbl.hasNumVec("speeds"))
    {
      s = tbl.getNumVec("speeds");
      for (int i=0; i<std::min<int>(3, s.size()); ++i)
        u[i] = s[i];
    }
    else
      throw Lucee::Except("AdvectionEquation::readInput: Must provide speeds in each direction");
  }

  void
  AdvectionEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  void
  AdvectionEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  void
  AdvectionEquation::flux(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
    double v[3];
    c.rotateVecToLocal(u, v);
    f[0] = v[0]*q[0];
  }

  void
  AdvectionEquation::speeds(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, double s[2])
  {
    double v[3];
    c.rotateVecToLocal(u, v);
    s[0] = v[0];
    s[1] = v[0];
  }

  void
  AdvectionEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double v[3];
    c.rotateVecToLocal(u, v);
    waves(0,0) = jump[0];
    s[0] = v[0];
  }

  void
  AdvectionEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v)
  {
    v[0] = q[0];
  }

  void
  AdvectionEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q)
  {
    q[0] = v[0];
  }
}

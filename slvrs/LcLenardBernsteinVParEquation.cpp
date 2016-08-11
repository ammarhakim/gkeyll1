/**
 * @file        LcLenardBernsteinVParEquation.cpp
 *
 * @brief	Lenard-Bernstein collision operator in v_parallel
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// Gkeyll includes
#include <LcLenardBernsteinVParEquation.h>

namespace Lucee
{
  const char *LenardBernsteinVParEquation::id = "LenardBernsteinVPar";

  LenardBernsteinVParEquation::LenardBernsteinVParEquation()
    : Lucee::HyperEquation(2,0)
  {
  }

  void
  LenardBernsteinVParEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
  }

  void
  LenardBernsteinVParEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
  }

  void
  LenardBernsteinVParEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
  }

  void
  LenardBernsteinVParEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
  }
  
  double
  LenardBernsteinVParEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    return 0.0;
  }
}

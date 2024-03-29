/**
 * @file	LcAdvectionEquation.cpp
 *
 * @brief	Advection equation
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

    fluxType = UPWIND_FLUX;
    if (tbl.hasString("fluxType"))
    {
      if (tbl.getString("fluxType") == "central")
        fluxType = CENTRAL_FLUX;
      else if (tbl.getString("fluxType") == "upwind")
        fluxType = UPWIND_FLUX;
      else if (tbl.getString("fluxType") == "downwind")
        fluxType = DOWNWIND_FLUX;
      else
      {
        Lucee::Except lce("AdvectionEquation::readInput: 'flux' must be one of ");
        lce << " 'central', 'upwind', or 'downwind'. Provided '" << tbl.getString("flux")
            << "' instead";
        throw lce;
      }
    }
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
  AdvectionEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    double v[3];
    c.rotateVecToLocal(u, v);
    f[0] = v[0]*q[0];
  }

  void
  AdvectionEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    double v[3];
    c.rotateVecToLocal(u, v);
    s[0] = v[0];
    s[1] = v[0];
  }

  void
  AdvectionEquation::primitive(const double* q, double* v) const
  {
    v[0] = q[0];
  }

  void
  AdvectionEquation::conserved(const double* v, double* q) const
  {
    q[0] = v[0];
  }

  void
  AdvectionEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
        const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,    
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double v[3];
    c.rotateVecToLocal(u, v);
    waves(0,0) = jump[0];
    s[0] = v[0];
  }

  double
  AdvectionEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    double v[3];
    c.rotateVecToLocal(u, v);

    if (fluxType == UPWIND_FLUX)
    {
      if (v[0] > 0)
        f[0] = v[0]*ql[0];
      else
        f[0] = v[0]*qr[0];
    }
    else if (fluxType == DOWNWIND_FLUX)
    {
      if (v[0] > 0)
        f[0] = v[0]*qr[0];
      else
        f[0] = v[0]*ql[0];
    }
    else
    {
// central flux
      f[0] = 0.5*v[0]*(ql[0]+qr[0]);
    }

    return v[0];
  }

  void
  AdvectionEquation::projectOnLeftEigenvectors(const Lucee::RectCoordSys& c,
    const double* q, const double* vec, double* coeff)
  { // left eigenmatrix is unit matrix
    coeff[0] = vec[0];
  }

  void
  AdvectionEquation::reconWithRightEigenvectors(const Lucee::RectCoordSys& c,
    const double* q, const double* coeff, double* vec)
  { // right eigenmatrix is unit matrix
    vec[0] = coeff[0];
  }
}

/**
 * @file	LcAdvectionEquationTemplated.cpp
 *
 * @brief	Advection equation for testing 4D basis functions.
 * Aligned rectangular coordinate system has been assumed.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAdvectionEquationTemplated.h>
#include <LcMathLib.h>

// std includes
#include <algorithm>

namespace Lucee
{
  template <> const char *AdvectionEquationTemplated<1>::id = "Advection1D";
  template <> const char *AdvectionEquationTemplated<2>::id = "Advection2D";
  template <> const char *AdvectionEquationTemplated<3>::id = "Advection3D";
  template <> const char *AdvectionEquationTemplated<4>::id = "Advection4D";
  template <> const char *AdvectionEquationTemplated<5>::id = "Advection5D";

  template <unsigned NDIM>
  AdvectionEquationTemplated<NDIM>::AdvectionEquationTemplated()
    : Lucee::HyperEquation(1, 1)
  {
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
    std::vector<double> s;
    
    // by default set speeds to 0.0
    for (int i=0; i<NDIM; i++)
      u[i] = 0.0;

    if (tbl.hasNumVec("speeds"))
    {
      s = tbl.getNumVec("speeds");
      for (int i=0; i<std::min<int>(NDIM, s.size()); ++i)
        u[i] = s[i];
    }
    else
      throw Lucee::Except("AdvectionEquationTemplated::readInput: Must provide speeds in each direction");

    fluxType = UPWIND_FLUX;
    if (tbl.hasString("fluxType"))
    {
      if (tbl.getString("fluxType") == "central")
        fluxType = CENTRAL_FLUX;
      else if (tbl.getString("fluxType") == "upwind")
        fluxType = UPWIND_FLUX;
      else
      {
        Lucee::Except lce("AdvectionEquationTemplated::readInput: 'flux' must be one of ");
        lce << " 'central' or 'upwind'. Provided '" << tbl.getString("flux")
            << "' instead";
        throw lce;
      }
    }

    lowerBound = 0.0; // by default ensure positive solutions
    if (tbl.hasNumber("lowerBound"))
      lowerBound = tbl.getNumber("lowerBound");
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    double v = u[c.getAlignmentDirection()];
    f[0] = v*q[0];
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    double v = u[c.getAlignmentDirection()];
    s[0] = v;
    s[1] = v;
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::primitive(const double* q, double* v) const
  {
    v[0] = q[0];
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::conserved(const double* v, double* q) const
  {
    q[0] = v[0];
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double v = u[c.getAlignmentDirection()];
    waves(0,0) = jump[0];
    s[0] = v;
  }

  template <unsigned NDIM>
  double
  AdvectionEquationTemplated<NDIM>::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    double v = u[c.getAlignmentDirection()];

    if (fluxType == UPWIND_FLUX)
    {
      if (v > 0)
        f[0] = v*ql[0];
      else
        f[0] = v*qr[0];
    }
    else
    {
// central flux
      f[0] = 0.5*v*(ql[0]+qr[0]);
    }

    return v;
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::projectOnLeftEigenvectors(const Lucee::RectCoordSys& c,
    const double* q, const double* vec, double* coeff)
  { // left eigenmatrix is unit matrix
    coeff[0] = vec[0];
  }

  template <unsigned NDIM>
  void
  AdvectionEquationTemplated<NDIM>::reconWithRightEigenvectors(const Lucee::RectCoordSys& c,
    const double* q, const double* coeff, double* vec)
  { // right eigenmatrix is unit matrix
    vec[0] = coeff[0];
  }

  template <unsigned NDIM>
  bool
  AdvectionEquationTemplated<NDIM>::isInvariantDomain(const double* q) const
  {
    if (Lucee::isNan(q[0]) || q[0] <= lowerBound)
      return false;
    return true;
  }

  // instantiations
  template class AdvectionEquationTemplated<1>;
  template class AdvectionEquationTemplated<2>;
  template class AdvectionEquationTemplated<3>;
  template class AdvectionEquationTemplated<4>;
  template class AdvectionEquationTemplated<5>;
}

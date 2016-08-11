/**
 * @file        LcDivEquation.cpp
 *
 * @brief	Auxilary equations with unit flux.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDivEquation.h>

// std includes
#include <iostream>

namespace Lucee
{
// set id for creators
  template <> const char *DivEquation<1>::id = "DivAuxFlux1D";
  template <> const char *DivEquation<2>::id = "DivAuxFlux2D";
  template <> const char *DivEquation<3>::id = "DivAuxFlux3D";

  template <unsigned NDIM>
  DivEquation<NDIM>::DivEquation()
    : Lucee::HyperEquation(1, 0)
  {
  }

  template <unsigned NDIM>
  void
  DivEquation<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);

    fluxType = FIVE_POINT;
// flux type
    if (tbl.hasString("fluxType"))
    {
      if (tbl.getString("fluxType") == "3-point")
        fluxType = THREE_POINT;
      else if (tbl.getString("fluxType") == "5-point")
        fluxType = FIVE_POINT;
      else
      {
        Lucee::Except lce("DivEquation::readInput: 'fluxType' must be one of ");
        lce << "'3-point' or '5-point'. '" << tbl.getString("fluxType")
            << "' specified instead";
        throw lce;
      }
    }
  }

  template <unsigned NDIM>
  void
  DivEquation<NDIM>::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  template <unsigned NDIM>
  void
  DivEquation<NDIM>::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  template <unsigned NDIM>
  void
  DivEquation<NDIM>::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    double wvec[3] = {0, 0, 0}, wvecNorm[3];
    for (unsigned i=0; i<NDIM; ++i)
      wvec[i] = auxVars[0][i];

// rotate to bring into normal-tangent space
    c.rotateVecToLocal(wvec, wvecNorm);
    f[0] = wvecNorm[0];
  }

  template <unsigned NDIM>
  double
  DivEquation<NDIM>::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    double wvecl[3] = {0, 0, 0}, wvecr[3] = {0, 0, 0};
    double wvecNorml[3], wvecNormr[3];
    for (unsigned i=0; i<NDIM; ++i)
    {
      wvecl[i] = auxVarsl[0][i];
      wvecr[i] = auxVarsr[0][i];
    }
// rotate to bring into normal-tangent space
    c.rotateVecToLocal(wvecl, wvecNorml);
    c.rotateVecToLocal(wvecr, wvecNormr);

    if (fluxType == THREE_POINT)
      f[0] = wvecNormr[0];
    else
      f[0] = 0.5*(wvecNorml[0] + wvecNormr[0]);

    return 0.0;
  }

// instantiations
  template class DivEquation<1>;
  template class DivEquation<2>;
  template class DivEquation<3>;
}

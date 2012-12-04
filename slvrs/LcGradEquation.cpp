/**
 * @file        LcGradEquation.cpp
 *
 * @brief	Auxilary equations with unit flux.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGradEquation.h>

// std includes
#include <iostream>

namespace Lucee
{
// set id for creators
  template <> const char *GradEquation<1>::id = "GradAuxFlux1D";
  template <> const char *GradEquation<2>::id = "GradAuxFlux2D";
  template <> const char *GradEquation<3>::id = "GradAuxFlux3D";


  template <unsigned NDIM>
  GradEquation<NDIM>::GradEquation()
    : Lucee::HyperEquation(NDIM, 0)
  {
  }

  template <unsigned NDIM>
  void
  GradEquation<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
// coefficient
    alpha = tbl.getNumber("coefficient");
  }

  template <unsigned NDIM>
  void
  GradEquation<NDIM>::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    double in3[3] = {0,0,0}, out3[3] = {0,0,0};
    for (unsigned i=0; i<NDIM; ++i)
      in3[i] = inQ[i];
    c.rotateVecToLocal(in3, out3);
    for (unsigned i=0; i<NDIM; ++i)
      outQ[i] = out3[i];
  }

  template <unsigned NDIM>
  void
  GradEquation<NDIM>::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    double in3[3] = {0,0,0}, out3[3] = {0,0,0};
    for (unsigned i=0; i<NDIM; ++i)
      in3[i] = inQ[i];
    c.rotateVecToGlobal(&in3[0], &out3[0]);
    for (unsigned i=0; i<NDIM; ++i)
      outQ[i] = out3[i];
  }

  template <unsigned NDIM>
  void
  GradEquation<NDIM>::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    f[0] = alpha*auxVars[0][0];
    for (unsigned i=1; i<NDIM; ++i)
      f[i] = 0.0;
  }

  template <unsigned NDIM>
  double
  GradEquation<NDIM>::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    f[0] = 0.5*alpha*(auxVarsl[0][0] + auxVarsr[0][0]);
    for (unsigned i=1; i<NDIM; ++i)
      f[i] = 0.0;
    return 0.0;
  }

// instantiations
  template class GradEquation<1>;
  template class GradEquation<2>;
  template class GradEquation<3>;
}

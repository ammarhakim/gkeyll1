/**
 * @file        LcAuxAdvectionEquation.cpp
 *
 * @brief	Auxilary equations with unit flux.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAuxAdvectionEquation.h>

// std includes
#include <iostream>

namespace Lucee
{
// set id for creators
  template <> const char *AuxAdvectionEquation<1>::id = "AuxAdvection1D";
  template <> const char *AuxAdvectionEquation<2>::id = "AuxAdvection2D";
  template <> const char *AuxAdvectionEquation<3>::id = "AuxAdvection3D";
  template <> const char *AuxAdvectionEquation<4>::id = "AuxAdvection4D";
  template <> const char *AuxAdvectionEquation<5>::id = "AuxAdvection5D";
  template <> const char *AuxAdvectionEquation<6>::id = "AuxAdvection6D";

  static void rotateVecToLocal(unsigned dir, const double inVec[4], double outVec[4])
  {
    if (dir == 0)
    {
      outVec[0] = inVec[0];
      outVec[1] = inVec[1];
      outVec[2] = inVec[2];
      outVec[3] = inVec[3];
    }
    else if (dir == 1)
    {
      outVec[0] = inVec[1];
      outVec[1] = -inVec[0];
      outVec[2] = inVec[2];
      outVec[3] = inVec[3];
    }
    else if (dir == 2)
    {
      outVec[0] = inVec[2];
      outVec[1] = inVec[1];
      outVec[2] = -inVec[0];
      outVec[3] = inVec[3];
    }
    else
    {
      outVec[0] = inVec[3];
      outVec[1] = inVec[1];
      outVec[2] = inVec[2];
      outVec[3] = -inVec[0];
    }
  }

  static void rotateVecToGlobal(unsigned dir, const double inVec[4], double outVec[4]) 
  {
    if (dir == 0)
    {
      outVec[0] = inVec[0];
      outVec[1] = inVec[1];
      outVec[2] = inVec[2];
      outVec[3] = inVec[3];
    }
    else if (dir == 1)
    {
      outVec[0] = -inVec[1];
      outVec[1] = inVec[0];
      outVec[2] = inVec[2];
      outVec[3] = inVec[3];
    }
    else if (dir == 2)
    {
      outVec[0] = -inVec[2];
      outVec[1] = inVec[1];
      outVec[2] = inVec[0];
      outVec[3] = inVec[3];
    }
    else
    {
      outVec[0] = -inVec[3];
      outVec[1] = inVec[1];
      outVec[2] = inVec[2];
      outVec[3] = inVec[0];
    }
  }

  template <unsigned NDIM>
  AuxAdvectionEquation<NDIM>::AuxAdvectionEquation()
    : Lucee::HyperEquation(1, 0)
  {
  }

  template <unsigned NDIM>
  void
  AuxAdvectionEquation<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  AuxAdvectionEquation<NDIM>::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  template <unsigned NDIM>
  void
  AuxAdvectionEquation<NDIM>::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0];
  }

  template <unsigned NDIM>
  void
  AuxAdvectionEquation<NDIM>::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    double uvec[6] = {0, 0, 0, 0, 0, 0}, uvecNorm[6];
    for (unsigned i=0; i<NDIM; ++i)
      uvec[i] = auxVars[0][i];

    unsigned dir = c.getAlignmentDirection();
// rotate to bring into normal-tangent space
    rotateVecToLocal(dir, uvec, uvecNorm);
// compute flux
    f[0] = q[0]*uvecNorm[0];
  }

  template <unsigned NDIM>
  double
  AuxAdvectionEquation<NDIM>::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    double uvecl[6] = {0, 0, 0, 0, 0, 0}, uvecr[6] = {0, 0, 0, 0, 0, 0};
    double uvecNorml[4], uvecNormr[4];
    for (unsigned i=0; i<NDIM; ++i)
    {
      uvecl[i] = auxVarsl[0][i];
      uvecr[i] = auxVarsr[0][i];
    }

    unsigned dir = c.getAlignmentDirection();
// rotate to bring into normal-tangent space
    rotateVecToLocal(dir, uvecl, uvecNorml);
    rotateVecToLocal(dir, uvecr, uvecNormr);

    double absMaxs = std::max(std::fabs(uvecNorml[0]), std::fabs(uvecNormr[0]));
// use Lax-fluxes
    f[0] = 0.5*(ql[0]*uvecNorml[0] + qr[0]*uvecNormr[0]) - 0.5*absMaxs*(qr[0]-ql[0]);

    return absMaxs;
  }

// instantiations
  template class AuxAdvectionEquation<1>;
  template class AuxAdvectionEquation<2>;
  template class AuxAdvectionEquation<3>;
  template class AuxAdvectionEquation<4>;
  template class AuxAdvectionEquation<5>;
  template class AuxAdvectionEquation<6>;
}

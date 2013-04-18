/**
 * @file	LcTenMomentFluidSource.cpp
 *
 * @brief       Compute source terms in pressure tensor equations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTenMomentFluidSource.h>

namespace Lucee
{
// set id for creators
  const char *TenMomentFluidSource::id = "TenMomentFluid";

// makes indexing easier
  static const unsigned RHO = 0;
  static const unsigned U1 = 1;
  static const unsigned U2 = 2;
  static const unsigned U3 = 3;
  static const unsigned P11 = 4;
  static const unsigned P12 = 5;
  static const unsigned P13 = 6;
  static const unsigned P22 = 7;
  static const unsigned P23 = 8;
  static const unsigned P33 = 9;

  TenMomentFluidSource::TenMomentFluidSource()
    : Lucee::PointSourceIfc(16, 10)
  { 
// takes in [rho, rho*u, rho*v, rho*w, P11, P12, P13, P22, P23 P33,
//   Ex, Ey, Ez, Bx, By, Bz] and computes
// sources for [rho, rho*u, rho*v, rho*w, P11, P12, P13, P22, P23, P33]
  }

  void
  TenMomentFluidSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);

    double charge = tbl.getNumber("charge");
    double mass = tbl.getNumber("mass");
    qbym = charge/mass;

    hasCollisions = false;
    if (tbl.hasNumber("collisionFrequency"))
    {
      hasCollisions = true;
      nu = tbl.getNumber("collisionFrequency");
    }
  }

  void
  TenMomentFluidSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
// takes in [rho, rho*u, rho*v, rho*w, P11, P12, P13, P22, P23 P33,
//   Ex, Ey, Ez, Bx, By, Bz]
    double r = this->getData(0);
    double u = this->getData(1)/r;
    double v = this->getData(2)/r;
    double w = this->getData(3)/r;
    double pxx = this->getData(4);
    double pxy = this->getData(5);
    double pxz = this->getData(6);
    double pyy = this->getData(7);
    double pyz = this->getData(8);
    double pzz = this->getData(9);
    double ex = this->getData(10);
    double ey = this->getData(11);
    double ez = this->getData(12);
    double bx = this->getData(13);
    double by = this->getData(14);
    double bz = this->getData(15);

    double re = qbym;

    src[RHO] = 0.0;
// momentum source terms
    src[U1] = r*re*(ex+v*bz-w*by);
    src[U2] = r*re*(ey+w*bx-u*bz);
    src[U3] = r*re*(ez+u*by-v*bx);
// pressure tensor source terms
    src[P11] = 2.0*r*re*u*ex+2.0*re*(bz*pxy-by*pxz);
    src[P12] = r*re*(u*ey+v*ex)+re*(bz*pyy-by*pyz-bz*pxx+bx*pxz);
    src[P13] = r*re*(u*ez+w*ex)+re*(bz*pyz+by*pxx-by*pzz-bx*pxy);
    src[P22] = 2.0*r*re*v*ey+2.0*re*(bx*pyz-bz*pxy);
    src[P23] = r*re*(v*ez+w*ey)+re*(by*pxy-bz*pxz+bx*pzz-bx*pyy);
    src[P33] = 2.0*r*re*w*ez+2.0*re*(by*pxz-bx*pyz);

    if (hasCollisions)
    {
// compute pressure tensor in fluid frame
      double fpxx = pxx-r*u*u;
      double fpxy = pxy-r*u*v;
      double fpxz = pxz-r*u*w;
      double fpyy = pyy-r*v*v;
      double fpyz = pyz-r*v*w;
      double fpzz = pzz-r*w*w;

// compute scalar pressure
      double pr = (fpxx+fpyy+fpzz)/3.0;

// add in collision terms
      src[P11] += nu*(pr-fpxx);
      src[P12] += nu*(0-fpxy);
      src[P13] += nu*(0-fpxz);
      src[P22] += nu*(pr-fpyy);
      src[P23] += nu*(0-fpyz);
      src[P33] += nu*(pr-fpzz);
    }
  }

  void
  TenMomentFluidSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
    throw Lucee::Except("TenMomentFluidSource::getSourceJac: Not implemented!");
  }
}

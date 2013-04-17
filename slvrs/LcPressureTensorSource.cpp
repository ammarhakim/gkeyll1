/**
 * @file	LcPressureTensorSource.cpp
 *
 * @brief       Compute source terms in pressure tensor equations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPressureTensorSource.h>

namespace Lucee
{
// set id for creators
  const char *PressureTensorSource::id = "PressureTensor";

  PressureTensorSource::PressureTensorSource()
    : Lucee::PointSourceIfc(16, 6)
  { 
// takes in [rho, rho*u, rho*v, rho*w, P11, P12, P13, P22, P23 P33,
//   Ex, Ey, Ez, Bx, By, Bz] and computes
// sources for [P11, P12, P13, P22, P23, P33]
  }

  void
  PressureTensorSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
// get charge and mass of species
    double charge = tbl.getNumber("charge");
    double mass = tbl.getNumber("mass");
    qbym = charge/mass;
  }

  void
  PressureTensorSource::getSource(double tm, const double loc[3], std::vector<double>& src)
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

// computes sources for [P11, P12, P13, P22, P23 P33]
    src[0] = 2.0*r*re*u*ex+2.0*re*(bz*pxy-by*pxz);
    src[1] = r*re*(u*ey+v*ex)+re*(bz*pyy-by*pyz-bz*pxx+bx*pxz);
    src[2] = r*re*(u*ez+w*ex)+re*(bz*pyz+by*pxx-by*pzz-bx*pxy);
    src[3] = 2.0*r*re*v*ey+2.0*re*(bx*pyz-bz*pxy);
    src[4] = r*re*(v*ez+w*ey)+re*(by*pxy-bz*pxz+bx*pzz-bx*pyy);
    src[5] = 2.0*r*re*w*ez+2.0*re*(by*pxz-bx*pyz);
  }

  void
  PressureTensorSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
    throw Lucee::Except("PressureTensorSource::getSourceJac: Not implemented!");
  }
}

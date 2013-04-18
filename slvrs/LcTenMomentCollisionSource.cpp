/**
 * @file	LcTenMomentCollisionSource.cpp
 *
 * @brief       Compute collisional source terms in 10-moment fluid equations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTenMomentCollisionSource.h>

namespace Lucee
{
// set id for creators
  const char *TenMomentCollisionSource::id = "TenMomentCollision";

  TenMomentCollisionSource::TenMomentCollisionSource()
    : Lucee::PointSourceIfc(10, 6)
  { 
// takes in [rho, rho*u, rho*v, rho*w, P11, P12, P13, P22, P23 P33] and computes
// sources for [P11, P12, P13, P22, P23, P33]
  }

  void
  TenMomentCollisionSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
    nu = tbl.getNumber("collisionFrequency");
  }

  void
  TenMomentCollisionSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
// takes in [rho, rho*u, rho*v, rho*w, P11, P12, P13, P22, P23 P33]
    double r = this->getData(0);
    double u = this->getData(1)/r;
    double v = this->getData(2)/r;
    double w = this->getData(3)/r;
    double pxx = this->getData(4)-r*u*u;
    double pxy = this->getData(5)-r*u*v;
    double pxz = this->getData(6)-r*u*w;
    double pyy = this->getData(7)-r*v*v;
    double pyz = this->getData(8)-r*v*w;
    double pzz = this->getData(9)-r*w*w;

// compute scalar pressure
    double pr = (pxx+pyy+pzz)/3.0;

// add in collision terms
    src[0] = nu*(pr-pxx);
    src[1] = nu*(0-pxy);
    src[2] = nu*(0-pxz);
    src[3] = nu*(pr-pyy);
    src[4] = nu*(0-pyz);
    src[5] = nu*(pr-pzz);
  }

  void
  TenMomentCollisionSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
    throw Lucee::Except("TenMomentCollisionSource::getSourceJac: Not implemented!");
  }
}

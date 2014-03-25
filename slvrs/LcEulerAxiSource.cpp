/**
 * @file	LcEulerAxiSource.cpp
 *
 * @brief	Source for computing Lorentz force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEulerAxiSource.h>

namespace Lucee
{
// set id for creators
  const char *EulerAxiSource::id = "EulerAxisymmetric";

  EulerAxiSource::EulerAxiSource()
    : Lucee::PointSourceIfc(5, 5)
  { 
// takes in [rho, rho*u, rho*v, rho*w, Er] and computes
// sources for [rho, rho*u, rho*v, rho*w, Er]
  }

  void
  EulerAxiSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
// read out gas gamma
    gasGamma = tbl.getNumber("gasGamma");
  }

  void
  EulerAxiSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
// takes in [rho, rho*u, rho*v, rho*w, Er]
    double rho = this->getData(0);
    double u = this->getData(1)/rho;
    double v = this->getData(2)/rho;
    double w = this->getData(3)/rho;
    double er = this->getData(4);

    double r = loc[0]; // radial coordinate
    double pr = (gasGamma-1)*(er - 0.5*rho*(u*u+v*v+w*w));

// computes sources for [rho, rho*u, rho*v, rho*w, Er]. These follow
// from expanding out the divergence of a vector and tensor in
// (r,theta,z) coordinates.
    src[0] = -rho*u/r;
    src[1] = -rho*u*u/r + rho*v*v/r;
    src[2] = -2*rho*u*v/r;
    src[3] = -rho*u*w/r;
    src[5] = -u*(er+pr)/r;
  }
}

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
  }

  void
  EulerAxiSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
// takes in [rho, rho*u, rho*v, rho*w, Er]
    double rho = this->getData(0);
    double rhou = this->getData(1);
    double rhov = this->getData(2);
    double rhow = this->getData(3);
    double er = this->getData(4);

    double r = loc[0]; // radial coordinate
// computes sources for [rho, rho*u, rho*v, rho*w, Er]
    src[0] = 0.0;
    src[1] = 0.0;
    src[2] = 0.0;
    src[3] = 0.0;
    src[5] = 0.0;
  }
}

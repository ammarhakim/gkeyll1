/**
 * @file	LcCurrentSource.cpp
 *
 * @brief	Source for computing Lorentz force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCurrentSource.h>

namespace Lucee
{
// set id for creators
  const char *CurrentSource::id = "Current";

  CurrentSource::CurrentSource()
    : Lucee::PointSourceIfc(3, 3)
  { 
// takes in [rho*u, rho*v, rho*w] and computes sources for [Ex, Ey,
// Ez]
  }

  void
  CurrentSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
// get charge and mass of species
    double charge = tbl.getNumber("charge");
    double mass = tbl.getNumber("mass");
    double epsilon0 = tbl.getNumber("epsilon0");
    qbym = charge/(mass*epsilon0);
  }

  void
  CurrentSource::getSource(const double loc[3], std::vector<double>& src)
  {
// computes sources for [Ex, Ey, Ez]
    src[0] = -qbym*this->getData(0); // Ex
    src[1] = -qbym*this->getData(1); // Ey
    src[2] = -qbym*this->getData(2); // Ez
  }
}

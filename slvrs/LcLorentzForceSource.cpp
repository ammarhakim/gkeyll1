/**
 * @file	LcLorentzForceSource.cpp
 *
 * @brief	Source for computing Lorentz force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLorentzForceSource.h>

namespace Lucee
{
// set id for creators
  const char *LorentzForceSource::id = "LorentzForce";

// indices into the source vector
  const static unsigned RHOU = 0;
  const static unsigned RHOV = 1;
  const static unsigned RHOW = 2;
  const static unsigned ER = 3;

  LorentzForceSource::LorentzForceSource()
    : Lucee::PointSourceIfc(10, 4)
  { 
// takes in [rho, rho*u, rho*v, rho*w, Ex, Ey, Ez, Bx, By, Bz] and computes
// sources for [rho*u, rho*v, rho*w, Er]
  }

  void
  LorentzForceSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
// get charge and mass of species
    double charge = tbl.getNumber("charge");
    double mass = tbl.getNumber("mass");
    qbym = charge/mass;
  }

  void
  LorentzForceSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
// takes in [rho, rho*u, rho*v, rho*w, Ex, Ey, Ez, Bx, By, Bz]
    double rho = this->getData(0);
    double rhou = this->getData(1);
    double rhov = this->getData(2);
    double rhow = this->getData(3);
    double ex = this->getData(4);
    double ey = this->getData(5);
    double ez = this->getData(6);
    double bx = this->getData(7);
    double by = this->getData(8);
    double bz = this->getData(9);

// computes sources for [rho*u, rho*v, rho*w, Er]
    src[0] = qbym*(rho*ex+rhov*bz-rhow*by); // x-momentum
    src[1] = qbym*(rho*ey+rhow*bx-rhou*bz); // y-momentum
    src[2] = qbym*(rho*ez+rhou*by-rhov*bx); // z-momentum
    src[3] = qbym*(rhou*ex+rhov*ey+rhow*ez); // energy
//     std::cout << "Lorentz source ..." << std::endl;
//     for (unsigned i=0; i<4; ++i)
//       std::cout << src[i] << std::endl;
//     std::cout << ".. done " << std::endl;
  }

  void
  LorentzForceSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
// takes in [rho, rho*u, rho*v, rho*w, Ex, Ey, Ez, Bx, By, Bz]
    double rho = this->getData(0);
    double rhou = this->getData(1);
    double rhov = this->getData(2);
    double rhow = this->getData(3);
    double ex = this->getData(4);
    double ey = this->getData(5);
    double ez = this->getData(6);
    double bx = this->getData(7);
    double by = this->getData(8);
    double bz = this->getData(9);

// row RHOU
    jac(RHOU, RHOV) = qbym*bz;
    jac(RHOU, RHOW) = -qbym*by;
// row RHOV
    jac(RHOV, RHOU) = -qbym*bz;
    jac(RHOV, RHOW) = qbym*bx;
// rho RHOW
    jac(RHOW, RHOU) = qbym*by;
    jac(RHOW, RHOV) = -qbym*bx;
// rho ER
    jac(ER, RHOU) = qbym*ex;
    jac(ER, RHOV) = qbym*ey;
    jac(ER, RHOW) = qbym*ez;
  }
}

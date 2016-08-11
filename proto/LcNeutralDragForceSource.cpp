/**
 * @file	LcNeutralDragForceSource.cpp
 *
 * @brief	Source for computing NeutralDrag force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNeutralDragForceSource.h>

namespace Lucee
{
// set id for creators
  const char *NeutralDragForceSource::id = "NeutralDragForce";

  NeutralDragForceSource::NeutralDragForceSource()
    : Lucee::PointSourceIfc(4, 4)
  { 
    // takes in [rho, rho*u, rho*v, rho*w,] and computes
    // sources for [rho*u, rho*v, rho*w, Er]
  }

  void
  NeutralDragForceSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
    // get chargeIn and mass of species
    if (tbl.hasNumVec("velocityNeut"))
    {
      velocityNeut = tbl.getNumVec("velocityNeut");
      if (velocityNeut.size() != 3)
        throw Lucee::Except(
          "NeutralDragForce::readInput: 'velocityNeut' table must have size of 3");
    }
    else
      throw Lucee::Except("NeutralDragForce::readInput: 'velocityNeut' table must be specified");

    if (tbl.hasNumber("nu"))
      nu = tbl.getNumber("nu");
    else
      throw Lucee::Except("NeutralDragForce::readInput: Collision frequency 'nu' must be specified");
  }

  void
  NeutralDragForceSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    double momentumRelU; //relative x-momentum
    double momentumRelV; //relative x-momentum
    double momentumRelW; //relative x-momentum
    // takes in [rho, rho*u, rho*v, rho*w]
    double rho = this->getData(0);
    double rhou = this->getData(1);
    double rhov = this->getData(2);
    double rhow = this->getData(3);

    // computes sources for [rho*u, rho*v, rho*w, Er]
    momentumRelU = rhou-rho*velocityNeut[0];
    momentumRelV = rhov-rho*velocityNeut[1];
    momentumRelW = rhow-rho*velocityNeut[2];
    src[0] = -nu*momentumRelU; // x-momentum
    src[1] = -nu*momentumRelV; // y-momentum
    src[2] = -nu*momentumRelW; // z-momentum
    src[3] = -0.5*nu*(momentumRelU*momentumRelU+momentumRelV*momentumRelV+momentumRelW*momentumRelW)/rho; // energy
  }

  void
  NeutralDragForceSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
  }
}

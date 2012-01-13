/**
 * @file	LcMaxwellTm2DUpdater.cpp
 *
 * @brief	Solver for transverse-magnetic Maxwell equations in 2D.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMathPhysConstants.h>
#include <LcMaxwellTm2DUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{
// constants to make indexing easier
  static const int BX = 0, BY = 1, EZ = 2, PSI = 3;

// set ids for module system
  const char *MaxwellTm2DUpdater::id = "MaxwellTm2D";

  MaxwellTm2DUpdater::MaxwellTm2DUpdater() 
    : Lucee::UpdaterIfc() {
  }

  void
  MaxwellTm2DUpdater::readInput(Lucee::LuaTable& tbl) {
// call base class method
    UpdaterIfc::readInput(tbl);
// speed of light
    c0 = Lucee::SPEED_OF_LIGHT;
    if (tbl.hasNumber("c0"))
      c0 = tbl.getNumber("c0");

// correction speed factor
    gamma = 1.0;
    if (tbl.hasNumber("gamma"))
      gamma = tbl.getNumber("gamma");
  }

  Lucee::UpdaterStatus
  MaxwellTm2DUpdater::update(double t) {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    double dx = grid.getDx(0), dy = grid.getDx(1);
    double dt = t-this->getCurrTime();

    double maxDtx = 0.5*dx/c0, maxDty = 0.5*dy/c0;
    double maxDt = maxDtx < maxDty ? maxDtx : maxDty; // minimum dt
    if (dt > maxDt)
// dt is too large
      return Lucee::UpdaterStatus(false, maxDt);

// compute various ratios needed in calculations
    double dtdx = 0.5*dt/dx, gc2dtdx = 0.5*gamma*c0*c0*dt/dx,
      gdtdx = 0.5*gamma*dt/dx, g2c2dtdx2 = 0.5*gamma*gamma*c0*c0*dt*dt/(dx*dx),
      c2dtdx = 0.5*c0*c0*dt/dx, c2dtdx2 = 0.5*c0*c0*dt*dt/(dx*dx);
    double dtdy = 0.5*dt/dy, gc2dtdy = 0.5*gamma*c0*c0*dt/dy,
      gdtdy = 0.5*gamma*dt/dy, g2c2dtdy2 = 0.5*gamma*gamma*c0*c0*dt*dt/(dy*dy),
      c2dtdy = 0.5*c0*c0*dt/dy, c2dtdy2 = 0.5*c0*c0*dt*dt/(dy*dy);

// get input/output arrays
    const Lucee::Field<2, double>& q = this->getInp<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& qNew = this->getOut<Lucee::Field<2, double> >(0);

// loop, updating solution
    Lucee::Region<2, int> localBox = grid.getLocalRegion();
    for (int i=localBox.getLower(0); i<localBox.getUpper(0); ++i)
    {
      for (int j=localBox.getLower(1); j<localBox.getUpper(1); ++j)
      {
// update Bx
        qNew(i,j,BX) = q(i,j,BX)
          - dtdy*( q(i,j+1,EZ) - q(i,j-1,EZ) )
          - gdtdx*( q(i+1,j,PSI)-q(i-1,j,PSI) )
          + c2dtdy2*( q(i,j+1,BX)-2*q(i,j,BX)+q(i,j-1,BX) )
          + g2c2dtdx2* (q(i+1,j,BX)-2*q(i,j,BX)+q(i-1,j,BX) );
// update By
        qNew(i,j,BY) = q(i,j,BY) 
          + dtdx*( q(i+1,j,EZ) - q(i-1,j,EZ) )
          - gdtdy*( q(i,j+1,PSI)-q(i,j-1,PSI) )
          + c2dtdx2*( q(i+1,j,BY)-2*q(i,j,BY)+q(i-1,j,BY) )
          + g2c2dtdy2* (q(i,j+1,BY)-2*q(i,j,BY)+q(i,j-1,BY) );
// update Ez
        qNew(i,j,EZ) = q(i,j,EZ)
          + c2dtdx*( q(i+1,j,BY)-q(i-1,j,BY) ) - c2dtdy*( q(i,j+1,BX)-q(i,j-1,BX) )
          + c2dtdx2*( q(i+1,j,EZ)-2*q(i,j,EZ)+q(i-1,j,EZ) )
          + c2dtdy2*( q(i,j+1,EZ)-2*q(i,j,EZ)+q(i,j-1,EZ) );
// update psi
//         qNew(i,j,PSI) = q(i,j,PSI)
//           - gc2dtdx*( q(i+1,j,BX)-q(i-1,j,BX) ) - gc2dtdy*( q(i,j+1,BY)-q(i,j-1,BY) )
//           + g2c2dtdx2*( q(i+1,j,PSI)-2*q(i,j,PSI)+q(i-1,j,PSI) )
//           + g2c2dtdy2*( q(i,j+1,PSI)-2*q(i,j,PSI)+q(i,j-1,PSI) );
        qNew(i,j,PSI) = q(i,j,PSI)
          - gc2dtdx*( q(i+1,j,BX)-q(i-1,j,BX) ) - gc2dtdy*( q(i,j+1,BY)-q(i,j-1,BY) );
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  MaxwellTm2DUpdater::declareTypes() {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}

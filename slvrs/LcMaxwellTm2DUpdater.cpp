/**
 * @file	LcMaxwellTm2DUpdater.cpp
 *
 * @brief	Solver for transverse-magnetic Maxwell equations in 2D.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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

namespace Lucee
{
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
// compute time-step

    return Lucee::UpdaterStatus();
  }

  void
  MaxwellTm2DUpdater::declareTypes() {
// expects 1 electric field and 1 magnetic field as input
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// expects 1 electric field and 1 magnetic field as output
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}

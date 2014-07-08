/**
 * @file	LcStairSteppedBcUpdater.cpp
 *
 * @brief	Base class for boundary conditions for stair-stepped boundaries.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcStairSteppedBcUpdater.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *StairSteppedBcUpdater<1>::id = "StairSteppedBc1D";
  template <> const char *StairSteppedBcUpdater<2>::id = "StairSteppedBc2D";
  template <> const char *StairSteppedBcUpdater<3>::id = "StairSteppedBc3D";

  template <unsigned NDIM>
  StairSteppedBcUpdater<NDIM>::StairSteppedBcUpdater() 
    : Lucee::UpdaterIfc() 
  {
  }

  template <unsigned NDIM>
  StairSteppedBcUpdater<NDIM>::~StairSteppedBcUpdater() 
  {
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl) 
  {
    UpdaterIfc::readInput(tbl);
// read direction to apply
    dir = tbl.getNumber("dir");

// get list of boundary conditions to apply
    if ( !tbl.hasTable("boundaryConditions") )
      Lucee::Except lce("StairSteppedBcUpdater::readInput: Must specify boundary conditions to apply!");

    Lucee::LuaTable bcTbl = tbl.getTable("boundaryConditions");
    bcList = bcTbl.template getAllObjects<Lucee::BoundaryCondition>();
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::initialize() 
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  StairSteppedBcUpdater<NDIM>::update(double t) 
  {
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::declareTypes() 
  {
// any number of output fields
    this->setLastInpVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class StairSteppedBcUpdater<1>;
  template class StairSteppedBcUpdater<2>;
  template class StairSteppedBcUpdater<3>;
}

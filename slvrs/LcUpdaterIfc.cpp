/**
 * @file	LcUpdaterIfc.cpp
 *
 * @brief	Base class for updaters in Lucee.
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
#include <LcExcept.h>
#include <LcSolverAssembly.h>
#include <LcUpdaterIfc.h>

// std includes
#include <limits>

namespace Lucee
{
// set module name
  const char *UpdaterIfc::id = "Updater";

  UpdaterIfc::UpdaterIfc()
    : Lucee::BasicObj(UpdaterIfc::id), grid(0)
  {
  }

  UpdaterIfc::~UpdaterIfc()
  {
  }

  void
  UpdaterIfc::readInput(Lucee::LuaTable& tbl)
  {
    throw Lucee::Except("UpdaterIfc::readInput: This method is not implemented");
  }

  void
  UpdaterIfc::initialize()
  {
  }

  void
  UpdaterIfc::setCurrTime(double tm) 
  {
    currTime = tm;
  }

  double
  UpdaterIfc::getCurrTime() const 
  {
    return currTime;
  }

  void
  UpdaterIfc::setGrid(const Lucee::GridIfc& grd)
  {
    grid = &grd;
  }

  void
  UpdaterIfc::setInpVars(const std::vector<const Lucee::DataStructIfc*>& dsl)
  {
    inpVars.resize(dsl.size());
     for (unsigned i=0; i<dsl.size(); ++i)
       inpVars[i] = dsl[i];
  }

  void
  UpdaterIfc::setOutVars(const std::vector<Lucee::DataStructIfc*>& dsl)
  {
    outVars.resize(dsl.size());
    for (unsigned i=0; i<dsl.size(); ++i)
      outVars[i] = dsl[i];
  }

  void
  UpdaterIfc::appendInpVarType(const std::type_info& type)
  {
    inpVarTypes.varTypes.push_back(&type);
  }

  void
  UpdaterIfc::appendOutVarType(const std::type_info& type)
  {
    outVarTypes.varTypes.push_back(&type);
  }

  void
  UpdaterIfc::setLastInpVarType(const std::type_info& type)
  {
    inpVarTypes.lastVarType = &type;
  }

  void
  UpdaterIfc::setLastOutVarType(const std::type_info& type)
  {
    outVarTypes.lastVarType = &type;
  }
}

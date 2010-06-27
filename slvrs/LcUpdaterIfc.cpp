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
#include <LcSolverAssembly.h>
#include <LcUpdaterIfc.h>

// std includes
#include <limits>

namespace Lucee
{
// set module name
  const char *UpdaterIfc::id = "Updater";

  UpdaterIfc::UpdaterIfc()
    : SolverIfc(UpdaterIfc::id)
  {
  }

  UpdaterIfc::~UpdaterIfc()
  {
  }

  void
  UpdaterIfc::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    SolverIfc::readInput(tbl);
  }

  void
  UpdaterIfc::buildData()
  {
  }

  void
  UpdaterIfc::buildAlgorithms()
  {
  }

  void
  UpdaterIfc::initialize()
  {
  }

  double
  UpdaterIfc::getSuggestedDt()
  {
    return std::numeric_limits<double>::max();
  }

  void
  UpdaterIfc::writeToFile(const std::string& baseName, unsigned d)
  { // updaters in general should not write anything to file
  }

  void
  UpdaterIfc::restoreFromFile(const std::string& baseName)
  { // updaters in general should not restore anything from file
  }

  void
  UpdaterIfc::finalize()
  {
  }

  void
  UpdaterIfc::setGrid(const Lucee::GridIfc& grd)
  {
    grid = &grd;
  }

  void
  UpdaterIfc::setInpVar(unsigned loc, const Lucee::DataStructIfc& ds)
  {
    inpVars[loc] = &ds;
  }

  void
  UpdaterIfc::setOutVar(unsigned loc, Lucee::DataStructIfc& ds)
  {
    outVars[loc] = &ds;
  }

  bool
  UpdaterIfc::typeCheck(const std::vector<std::string>& inp,
    const std::vector<std::string>& out) const
  {
    return true;
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

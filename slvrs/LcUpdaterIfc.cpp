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
  UpdaterIfc::setInpVarNames(const std::vector<std::string>& nms)
  {
  }

  void
  UpdaterIfc::setOutVarNames(const std::vector<std::string>& nms)
  {
  }
}

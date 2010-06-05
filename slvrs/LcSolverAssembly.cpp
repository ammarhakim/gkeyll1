/**
 * @file	LcSolverAssembly.cpp
 *
 * @brief	Solver that assembles updaters to create a simulation. *
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

namespace Lucee
{
  const char *SolverAssembly::id = "Assembly";
  
  SolverAssembly::SolverAssembly()
    : SolverIfc(SolverAssembly::id)
  {
  }

  SolverAssembly::~SolverAssembly()
  {
// delete all grids
    std::map<std::string, Lucee::GridBase*>::iterator gItr;
    for (gItr = gridMap.begin(); gItr != gridMap.end(); ++gItr)
      delete gItr->second;
    gridMap.clear();
  }
  
  void 
  SolverAssembly::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  SolverAssembly::buildData()
  {
  }

  void 
  SolverAssembly::buildAlgorithms()
  {
  }

  void
  SolverAssembly::initialize()
  {
  }
  
  void 
  SolverAssembly::writeToFile(const std::string& baseName, unsigned d)
  {
  }

  int
  SolverAssembly::advance(double t)
  {
    return 0;
  }

  void
  SolverAssembly::restoreFromFile(const std::string& baseName)
  {
  }

  void
  SolverAssembly::finalize()
  {
  }
}

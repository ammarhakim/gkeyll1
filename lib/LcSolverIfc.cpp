/**
 * @file	LcSolverIfc.cpp
 *
 * @brief	Interface class for Lucee solvers.
 *
 * @version	$Id: LcSolverIfc.cpp 312 2010-03-02 18:37:24Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSolverIfc.h>

// std includes
#include <string>

namespace Lucee
{
// set madule name
  const char *SolverIfc::id = "Solver";

  SolverIfc::SolverIfc(const std::string& nm)
    : Lucee::BasicObj(nm)
  {
  }

  SolverIfc::~SolverIfc()
  {
  }

  void
  SolverIfc::setCurrTime(double tm) 
  {
    currTime = tm;
  }

  void
  SolverIfc::setNumOutFrames(unsigned n)
  {
    numOutFrames = n;
  }

  void
  SolverIfc::readInput(Lucee::LuaTable& tbl)
  {
  }
}

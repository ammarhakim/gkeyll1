/**
 * @file	LcFunctionBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFunctionBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *FunctionBoundaryCondition::id = "Function";

  void
  FunctionBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);
  }

  void
  FunctionBoundaryCondition::applyBc(double tm, const double loc[3],
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
  }
}

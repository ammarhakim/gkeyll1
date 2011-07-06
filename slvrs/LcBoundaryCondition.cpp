/**
 * @file	LcBoundaryCondition.cpp
 *
 * @brief	Base class for applying boundary condition.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>

namespace Lucee
{
// set module name
  const char *BoundaryCondition::id = "BoundaryCondition";

  BoundaryCondition::BoundaryCondition()
    : Lucee::BasicObj(BoundaryCondition::id)
  {
  }

  void
  BoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
  }
}

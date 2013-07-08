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
    Lucee::BasicObj::readInput(tbl);
    if (tbl.hasNumVec("components"))
    {
// get list of components to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("components");
      for (unsigned i=0; i<cd.size(); ++i)
        components.push_back( (int) cd[i] );
    }
  }

  void
  BoundaryCondition::setComponents(const std::vector<unsigned>& cd)
  {
    components.clear();
    for (unsigned i=0; i<cd.size(); ++i)
      components.push_back( (int) cd[i] );
  }
}

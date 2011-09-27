/**
 * @file	LcConstBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcConstBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *ConstBoundaryCondition::id = "Const";

  void
  ConstBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);
// read out factors for multiplications
    if (tbl.hasNumVec("values"))
    {
      values = tbl.getNumVec("values");
      if (values.size() != components.size())
        throw Lucee::Except(
          "ConstBoundaryCondition::readInput: If 'values' table is specified it must have same size as 'components' table");
    }
    else
    {
      values.resize(components.size());
      for (unsigned i=0; i<values.size(); ++i) values[i] = 0.0;
    }
  }

  void
  ConstBoundaryCondition::applyBc(
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// just copy data over
    for (unsigned i=0; i<components.size(); ++i)
      qbc[components[i]] = values[i];
  }
}

/**
 * @file	LcCopyBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCopyBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *CopyBoundaryCondition::id = "Copy";

  void
  CopyBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);
// read out factors for multiplications
    if (tbl.hasNumVec("fact"))
    {
      fact = tbl.getNumVec("fact");
      if (fact.size() != components.size())
        throw Lucee::Except(
          "CopyBoundaryCondition::readInput: If 'fact' table is specified it must have same size as 'components' table");
    }
    else
    {
      fact.resize(components.size());
      for (unsigned i=0; i<fact.size(); ++i) fact[i] = 1.0;
    }
  }

  void
  CopyBoundaryCondition::applyBc(
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// just copy data over
    for (unsigned i=0; i<components.size(); ++i)
      qbc[components[i]] = fact[i]*qin[components[i]];
  }
}

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
      if (fact.size() != this->numComponents())
        throw Lucee::Except(
          "CopyBoundaryCondition::readInput: If 'fact' table is specified it must have same size as 'components' table");
    }
    else
    {
      fact.resize(this->numComponents());
      for (unsigned i=0; i<fact.size(); ++i) fact[i] = 1.0;
    }
  }

  void
  CopyBoundaryCondition::applyBc(double tm, const double loc[3],
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// just copy data over
    for (unsigned i=0; i<this->numComponents(); ++i)
      qbc[this->component(i)] = fact[i]*qin[this->component(i)];
  }
}

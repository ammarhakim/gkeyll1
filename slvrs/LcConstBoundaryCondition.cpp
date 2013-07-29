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
      if (values.size() != this->numComponents())
        throw Lucee::Except(
          "ConstBoundaryCondition::readInput: 'values' table must have same size as 'components' table");
    }
    else
    {
      throw Lucee::Except("ConstBoundaryCondition::readInput: 'values' table must be specified");
    }
  }

  void
  ConstBoundaryCondition::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// just copy data over
    for (unsigned i=0; i<this->numComponents(); ++i)
      qbc[this->component(i)] = values[i];
  }
}

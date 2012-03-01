/**
 * @file	LcNodalFiniteElementIfc.cpp
 *
 * @brief       Interface class for a reference nodal finite element.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGridIfc.h>
#include <LcNodalFiniteElementIfc.h>

namespace Lucee
{
// set module name
  const char *NodalFiniteElementIfc::id = "NodalFiniteElement";

  void
  NodalFiniteElementIfc::readInput(Lucee::LuaTable& tbl)
  {
// read in grid
    if (tbl.hasObject<Lucee::GridIfc>("onGrid"))
      grid = &tbl.getObjectAsBase<Lucee::GridIfc>("onGrid");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify grid using 'onGrid'");
  }

  double
  NodalFiniteElementIfc::evalBasis(unsigned n, double x, double y) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::evalBasis: Not implemented!");
    return 0;
  }

  void
  NodalFiniteElementIfc::getMassMatrix(Lucee::Matrix<double> NjNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getMassMatrix: Not implemented!");
  }

  void
  NodalFiniteElementIfc::getStiffnessMatrix(Lucee::Matrix<double> DNjDNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getStiffnessMatrix: Not implemented!");
  }

  NodalFiniteElementIfc::NodalFiniteElementIfc(unsigned numNodes)
    : numNodes(numNodes)
  {
  }
}

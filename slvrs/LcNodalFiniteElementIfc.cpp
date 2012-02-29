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
#include <LcNodalFiniteElementIfc.h>

namespace Lucee
{
// set module name
  const char *NodalFiniteElementIfc::id = "NodalFiniteElement";

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

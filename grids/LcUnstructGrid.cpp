/**
 * @file	LcUnstructGrid.cpp
 *
 * @brief	Unstructured grid class.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructGrid.h>

namespace Lucee
{
  template <typename REAL>
  UnstructGrid<REAL>::UnstructGrid()
    : ndim(3), ddprime(4*4), connectivity(4*4)
  {
// initialize ddprime as no connections have been made
    for (unsigned i=0; i<16; ++i)
      ddprime[i] = false;
  }

  template <typename REAL>
  void
  UnstructGrid<REAL>::constructFromCreator(const Lucee::UnstructGridCreator<REAL>& ctor)
  {
    ndim = ctor.getDim();
    ctor.fillWithGeometry(geometry);
    ctor.fillWithConnectivity(connectivity[4*ndim+0]);
    ddprime[4*ndim+0] = true; // set flag as now ndim->0 connectivity is stored
  }

// instantiations
  template class UnstructGrid<float>;
  template class UnstructGrid<double>;
}

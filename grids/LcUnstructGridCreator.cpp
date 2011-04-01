/**
 * @file	LcUnstructGridCreator.cpp
 *
 * @brief	Unstructured grid creator class.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructGridCreator.h>

namespace Lucee
{
  template <typename REAL>
  UnstructGridCreator<REAL>::UnstructGridCreator(unsigned ndim)
    : ndim(ndim)
  {
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setNumVertices(unsigned nv)
  {
    vc.reset(nv);
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setNumCells(unsigned nc)
  {
    c2v.reset(nc);
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::appendVertex(double xv[3])
  {
  }

// instantiations
  template class UnstructGridCreator<float>;
  template class UnstructGridCreator<double>;
}

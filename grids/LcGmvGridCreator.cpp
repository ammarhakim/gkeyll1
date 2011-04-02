/**
 * @file	LcGmvGridCreator.cpp
 *
 * @brief	Creator for grids defined in GMV format.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGmvGridCreator.h>

namespace Lucee
{
  template <typename REAL>
  GmvGridCreator<REAL>::GmvGridCreator(unsigned ndim, std::istream& gmvStrm)
    : UnstructGridCreator<REAL>(ndim)
  {
  }

// instantiations
  template class GmvGridCreator<float>;
  template class GmvGridCreator<double>;
}

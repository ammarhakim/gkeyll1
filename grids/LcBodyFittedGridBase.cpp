/**
 * @file	LcBodyFittedGridBase.cpp
 *
 * @brief	Base class for body fitted grid in arbitrary dimensions.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBodyFittedGridBase.h>

namespace Lucee
{
  template <unsigned NDIM>
  BodyFittedGridBase<NDIM>::BodyFittedGridBase(const Lucee::Region<NDIM, int>& localBox,
    const Lucee::Region<NDIM, int>& globalBox,
    const Lucee::Region<NDIM, double>& compSpace)
    : localBox(localBox), globalBox(globalBox), compSpace(compSpace)
  {
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  BodyFittedGridBase<NDIM>::getGlobalBox() const 
  { return globalBox; }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  BodyFittedGridBase<NDIM>::getLocalBox() const 
  { return localBox; }

  template <unsigned NDIM>
  Lucee::Region<NDIM, double>
  BodyFittedGridBase<NDIM>::getComputationalSpace() const 
  { return compSpace; }

// instantiations
  template class BodyFittedGridBase<1>;
  template class BodyFittedGridBase<2>;
  template class BodyFittedGridBase<3>;
}

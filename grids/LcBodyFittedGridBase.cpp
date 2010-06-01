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
  BodyFittedGridBase<NDIM>::~BodyFittedGridBase()
  {
  }

  template <unsigned NDIM>
  unsigned
  BodyFittedGridBase<NDIM>::getNumCells(unsigned dir) const
  {
    return globalBox.getShape(dir);
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  BodyFittedGridBase<NDIM>::getGlobalBox() const 
  { 
    return globalBox; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  BodyFittedGridBase<NDIM>::getLocalBox() const 
  { 
    return localBox; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, double>
  BodyFittedGridBase<NDIM>::getComputationalSpace() const 
  { 
    return compSpace; 
  }

  template <unsigned NDIM>
  void
  BodyFittedGridBase<NDIM>::setIndex(int i)
  {
    currIdx[0] = i;
  }
  
  template <unsigned NDIM>
  void
  BodyFittedGridBase<NDIM>::setIndex(int i, int j)
  {
    currIdx[0] = i;
    currIdx[1] = j;
  }

  template <unsigned NDIM>
  void
  BodyFittedGridBase<NDIM>::setIndex(int i, int j, int k)
  {
    currIdx[0] = i;
    currIdx[1] = j;
    currIdx[2] = k;
  }

  template <unsigned NDIM>
  void
  BodyFittedGridBase<NDIM>::setIndex(const int idx[NDIM])
  {
    for (unsigned i=0; i<NDIM; ++i)
      currIdx[i] = idx[i];
  }

// instantiations
  template class BodyFittedGridBase<1>;
  template class BodyFittedGridBase<2>;
  template class BodyFittedGridBase<3>;
}

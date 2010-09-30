/**
 * @file	LcStructuredGridBase.cpp
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
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <unsigned NDIM>
  StructuredGridBase<NDIM>::~StructuredGridBase()
  {
  }

  template <unsigned NDIM>
  unsigned
  StructuredGridBase<NDIM>::getNumCells(unsigned dir) const
  {
    return globalBox.getShape(dir);
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  StructuredGridBase<NDIM>::getGlobalBox() const 
  { 
    return globalBox; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  StructuredGridBase<NDIM>::getLocalBox() const 
  { 
    return localBox; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, double>
  StructuredGridBase<NDIM>::getComputationalSpace() const 
  { 
    return compSpace; 
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(int i) const
  {
    currIdx[0] = i;
  }
  
  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(int i, int j) const
  {
    currIdx[0] = i;
    currIdx[1] = j;
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(int i, int j, int k) const
  {
    currIdx[0] = i;
    currIdx[1] = j;
    currIdx[2] = k;
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(const int idx[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      currIdx[i] = idx[i];
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase()
    : localBox(&Lucee::FixedVector<NDIM, int>(1)[0]),
      globalBox(&Lucee::FixedVector<NDIM, int>(1)[0]),
      compSpace(&Lucee::FixedVector<NDIM, double>(1.0)[0])
  {
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase(const Lucee::Region<NDIM, int>& localBox,
    const Lucee::Region<NDIM, int>& globalBox,
    const Lucee::Region<NDIM, double>& compSpace)
    : localBox(localBox), globalBox(globalBox), compSpace(compSpace)
  {
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>&
  StructuredGridBase<NDIM>::operator=(const StructuredGridBase<NDIM>& sg)
  {
    if (&sg == this)
      return *this;

    for (unsigned i=0; i<NDIM; ++i)
      currIdx[i] = sg.currIdx[i];
    localBox = sg.localBox;
    globalBox = sg.globalBox;
    compSpace = sg.compSpace;
    
    return *this;
  }

// instantiations
  template class StructuredGridBase<1>;
  template class StructuredGridBase<2>;
  template class StructuredGridBase<3>;
}

/**
 * @file	LcRectCartGrid.cpp
 *
 * @brief	A rectangular cartesian grid.
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
#include <LcRectCartGrid.h>

namespace Lucee
{
  template <unsigned NDIM>
  RectCartGrid<NDIM>::RectCartGrid(const Lucee::Region<NDIM, int>& localBox,
    const Lucee::Region<NDIM, int>& globalBox,
    const Lucee::Region<NDIM, double>& physBox) 
    : Lucee::BodyFittedGridBase<NDIM>(localBox, globalBox, physBox)
  {
  }

// instantiations
  template class RectCartGrid<1>;
  template class RectCartGrid<2>;
  template class RectCartGrid<3>;
}

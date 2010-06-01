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
    for (unsigned i=0; i<3; ++i)
      dx[i] = 1.0;
    for (unsigned i=0; i<NDIM; ++i)
      physBox.getShape(i)/globalBox.getShape(i);
    cellVolume = 1.0;
    for (unsigned i=0; i<NDIM; ++i)
      cellVolume = cellVolume*dx[i];
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getCentriod(double xc[3]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      xc[i] = (currIdx[i]-globalBox.getLower(i) + 0.5)*dx[i];
    for (unsigned i=NDIM; i<3; ++i)
      xc[i] = 0.0;
  }

  template <unsigned NDIM>
  double
  RectCartGrid<NDIM>::getVolume() const
  {
    return cellVolume;
  }

  template <unsigned NDIM>
  double
  RectCartGrid<NDIM>::getSurfArea(unsigned dir) const
  {
    if (dir==0)
      return dx[1];
    else if (dir==1)
      return dx[0];
    else if (dir=2)
      retun dx[2];
    return 1.0;
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getSurfNormal(unsigned dir, double norm[3]) const
  {
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getSurfTangents(unsigned dir, double tan1[3], double tan2[3]) const
  {
  }

// instantiations
  template class RectCartGrid<1>;
  template class RectCartGrid<2>;
  template class RectCartGrid<3>;
}

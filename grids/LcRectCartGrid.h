/**
 * @file	LcRectCartGrid.h
 *
 * @brief	A rectangular cartesian grid.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_RECT_CART_GRID_H
#define LC_RECT_CART_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBodyFittedGridBase.h>

namespace Lucee
{
  template <unsigned NDIM>
  class RectCartGrid : public Lucee::BodyFittedGridBase<NDIM>
  {
    public:
/**
 * Create a new rectangular Cartesian grid on specified region. In
 * serial the local and global boxes coincide. In parallel, the
 * globalBox represents the full grid, while the localBox represents
 * the portion of the grid handled by the rank the grid lives on.
 *
 * To get data from the grid, first set the index into the grid by
 * using the setIndex() method. Then access the needed data for that
 * cell. For structured grids, the edges on the lower side (left,
 * bottom, back) are labeled by the cell index. Further, the lower
 * left corner of each cell is also labelled by the cell index.
 *
 * @param localBox Local index region for this grid.
 * @param globalBox Global index region for this grid.
 * @param physBox Rectangular region in space.
 */
      RectCartGrid(const Lucee::Region<NDIM, int>& localBox,
        const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& physBox);
  };
}

#endif // LC_RECT_CART_GRID_H

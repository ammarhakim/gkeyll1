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

/**
 * Return coordinates in physical space of cell centroid.
 *
 * @param xc On output, centroid of cell.
 */
      virtual void getCentriod(double xc[3]) const;

/**
 * Return volume of cell.
 *
 * @return Cell volume.
 */
      virtual double getVolume() const;

/**
 * Return physical surface area of face perpendicular (in
 * computational space) to specified direction.
 *
 * @param dir Direction perpendicular to face.
 * @return surface area of face.
 */
      virtual double getSurfArea(unsigned dir) const;

/**
 * Return unit normal to face perpendicular (in computational space)
 * to specified direction.
 *
 * @param dir Direction perpendicular to face.
 * @param norm On output, normal to face.
 */
      virtual void getSurfNormal(unsigned dir, double norm[3]) const;

/**
 * Return two mutually perpendicular unit tangents to face
 * perpendicular (in computational space) to specified direction. If
 * 'norm' is the surface normal then 'tan1 x tan2 = norm'.
 *
 * @param dir Direction perpendicular to face.
 * @param tan1 On output, first tangent to face.
 * @param tan2 On output, second tangent to face.
 */
      virtual void getSurfTangents(unsigned dir, double tan1[3], double tan2[3]) const;

    private:
/** Grid spacing in each direction */      
      double dx[3];
/** Volume of each cell */
      double cellVolume;
  };
}

#endif // LC_RECT_CART_GRID_H

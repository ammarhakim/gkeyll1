/**
 * @file	LcBodyFittedGridBase.h
 *
 * @brief	Base class for body fitted grid in arbitrary dimensions.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_BODY_FITTED_GRID_BASE_H
#define LC_BODY_FITTED_GRID_BASE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRegion.h>

namespace Lucee
{
/**
 * Body-fitted grid class in arbitrary dimensions. Provides an class
 * to represent a single-block rectangular body-fitted grid.
 */
  template <unsigned NDIM>
  class BodyFittedGridBase
  {
    public:
/**
 * Destroy object.
 */
      virtual ~BodyFittedGridBase();

/**
 * Get number of cells in given direction.
 *
 * @param dir Direction in which number of cells is needed.
 * @return Number of cells.
 */
      unsigned getNumCells(unsigned dir) const;

/**
 * Return global region indexed by grid.
 *
 * @return global region indexed by grid.
 */
      Lucee::Region<NDIM, int> getGlobalBox() const;

/**
 * Return local region indexed by grid.
 *
 * @return  local region indexed by grid.
 */
      Lucee::Region<NDIM, int> getLocalBox() const;

/**
 * Return region in computational space
 *
 * @return region in computational space
 */
      Lucee::Region<NDIM, double> getComputationalSpace() const;

/**
 * Set the current cell location in grid to (i).
 *
 * @param i Index location into grid.
 */
      void setIndex(int i);

/**
 * Set the current cell location in grid to (i, j).
 *
 * @param i Index location into grid.
 * @param j Index location into grid.
 */
      void setIndex(int i, int j);

/**
 * Set the current cell location in grid to (i, j, k).
 *
 * @param i Index location into grid.
 * @param j Index location into grid.
 * @param k Index location into grid.
 */
      void setIndex(int i, int j, int k);

/**
 * Set the current cell location in grid to specified index.
 *
 * @param idx Index location into grid.
 */
      void setIndex(const int idx[3]);

/**
 * Return coordinates in physical space of cell centroid.
 *
 * @param xc On output, centroid of cell.
 */
      virtual void getCentriod(double xc[3]) const = 0;

/**
 * Return volume of cell.
 *
 * @return Cell volume.
 */
      virtual double getVolume() const = 0;

/**
 * Return physical surface area of face perpendicular (in
 * computational space) to specified direction.
 *
 * @param dir Direction perpendicular to face.
 * @return surface area of face.
 */
      virtual double getSurfArea(unsigned dir) const = 0;

/**
 * Get the coordinate system attached to the surface perpendicular (in
 * computational space) to specified direction. The corrdinate system
 * is defined by three unit vectors, norm, tan1 and tan2. The norm
 * vector is normal to the surface and points into the cell. The
 * vectors tan1 and tan2 lie in the plane of the surface. The system
 * is such that tan1 x tan2 = norm.
 *
 * @param dir Direction perpendicular to face.
 * @param norm On output, normal to face.
 * @param tan1 On output, first tangent to face.
 * @param tan2 On output, second tangent to face.
 * 
 */
      virtual void getSurfCoordSys(unsigned dir, double norm[3],
        double tan1[3], double tan2[3]) const = 0;

    protected:
/**
 * Create a new body-fitted grid on specified region. In serial the
 * local and global boxes coincide. In parallel, the globalBox
 * represents the full grid, while the localBox represents the portion
 * of the grid handled by the rank the grid lives on.
 *
 * To get data from the grid, first set the index into the grid by
 * using the setIndex() method. Then access the needed data for that
 * cell. For structured grids, the edges on the lower side (left,
 * bottom, back) are labeled by the cell index. Further, the lower
 * left corner of each cell is also labelled by the cell index.
 *
 * @param localBox Local index region for this grid.
 * @param globalBox Global index region for this grid.
 * @param compSpace Region in computation space.
 */
      BodyFittedGridBase(const Lucee::Region<NDIM, int>& localBox,
        const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& compSpace);

/** Index into current cell */
      int currIdx[NDIM];
/** Local region indexed by grid */
      Lucee::Region<NDIM, int> localBox;
/** Global region indexed by grid */
      Lucee::Region<NDIM, int> globalBox;
/** Global region spanned by grid in computational space */
      Lucee::Region<NDIM, double> compSpace;
  };
}

#endif //  LC_BODY_FITTED_GRID_H

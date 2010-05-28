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
 * Create a new body-fitted grid on specified region. In serial the
 * local and global boxes coincide. In parallel, the globalBox
 * represents the full grid, while the localBox represents the portion
 * of the grid handled by the rank the grid lives on.
 *
 * @param localBox Local index region for this grid.
 * @param globalBox Global index region for this grid.
 * @param compSpace Region in computation space.
 */
      BodyFittedGridBase(const Lucee::Region<NDIM, int>& localBox,
        const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& compSpace);

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

    private:
/** Local region indexed by grid */
      Lucee::Region<NDIM, int> localBox;
/** Global region indexed by grid */
      Lucee::Region<NDIM, int> globalBox;
/** Region spanned by grid in computational space */
      Lucee::Region<NDIM, double> compSpace;
  };
}

#endif //  LC_BODY_FITTED_GRID_H

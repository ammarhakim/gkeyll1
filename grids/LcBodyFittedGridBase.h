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
 * Create a new body-fitted grid on specified region.
 *
 * @param globalBox Global index region for this grid.
 * @param compSpace Region in computation space.
 */
      BodyFittedGridBase(const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& compSpace);

    private:
/** Region indexed by grid */
      Lucee::Region<NDIM, int> globalBox;
/** Region spanned by grid in computational space */
      Lucee::Region<NDIM, double> compSpace;
  };
}

#endif //  LC_BODY_FITTED_GRID_H

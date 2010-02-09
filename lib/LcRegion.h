/**
 * @file	LcRegion.h
 *
 * @brief	Region in an N-dimensional space.
 *
 * @version	$Id: LcRegion.h 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_REGION_H
#define LC_REGION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFixedVector.h>

namespace Lucee
{
/**
 * A region in an N-dimension space represents a rectangular
 * hyper-box, closed on the lower end but open on the upper end. For
 * example, a two-dimensional region [a,b) X [c,d) has volume
 * (b-a)*(d-c).
 */
  template <unsigned NDIM, typename T>
  class Region
  {
    public:
/**
 * Create a region object with given shape. Lower limits are assumed
 * to be zero.
 *
 * @param shape Shape of region.
 */
      Region(T shape[NDIM]);

/**
 * Create a region object with given lower and upper bounds.
 *
 * @param lower Lower limits of region.
 * @param upper Upper limits of region.
 */
      Region(T lower[NDIM], T upper[NDIM]);

/**
 * Return lower bound of region in a given direction.
 *
 * @return Lower bound of region.
 */
      T getLower(unsigned i) const { return lower[i]; }

/**
 * Return upper bound of region in a given direction.
 *
 * @return Upper bound of region.
 */
      T getUpper(unsigned i) const { return upper[i]; }

/**
 * Return shape of region in a given direction.
 *
 * @return Shape of region.
 */
      T getShape(unsigned i) const { return upper[i]-lower[i]; }

/**
 * Get volume of the box.
 *
 * @return volume of box.
 */
      T getVolume() const { return volume; }

/**
 * Get lower bounds of box.
 *
 * @return lower bounds of box.
 */
      FixedVector<NDIM, T> getLower() const 
      { 
        return FixedVector<NDIM, T>(lower);
      }

/**
 * Get upper bounds of box.
 *
 * @return upper bounds of box.
 */
      FixedVector<NDIM, T> getUpper() const
      {
        return FixedVector<NDIM, T>(upper);
      }

/**
 * Check if point is inside region.
 *
 * @param point Point to check.
 */
      bool isInside(T point[NDIM]) const;

    private:
/** Lower and upper coordinates of region */
      T lower[NDIM], upper[NDIM];
/** Volume of region */
      T volume;
  };

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(T shape[NDIM])
  {
    volume = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      lower[i] = 0;
      upper[i] = shape[i];
      volume *= shape[i];
    }
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(T lo[NDIM], T up[NDIM])
  {
    volume = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      lower[i] = lo[i];
      upper[i] = up[i];
      volume *= (up[i]-lo[i]);
    }
  }

  template <unsigned NDIM, typename T>
  bool
  Region<NDIM, T>::isInside(T point[NDIM]) const
  {
    bool isIn = false;
    for (unsigned i=0; i<NDIM; ++i)
      if ((point[i] < lower[i]) || (point[i] >= upper[i]))
        return false;
    return true;
  }
}

#endif // LC_REGION_H

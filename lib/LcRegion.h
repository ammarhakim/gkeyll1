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
 * example, a two-dimensional region can be represented by the
 * Cartesian product [a,b) X [c,d) has volume (b-a)*(d-c).
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
 * Create a region object from an existing one.
 *
 * @param rgn Region object to copy.
 */
      Region(const Region<NDIM, T>& rgn);

/**
 * Copy a region object from an existing one.
 *
 * @param rgn Region object to copy.
 * @return reference to copy.
 */
      Region<NDIM, T>& operator=(const Region<NDIM, T>& rgn);

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
      FixedVector<NDIM, T> getLower() const { return lower; }

/**
 * Get upper bounds of box.
 *
 * @return upper bounds of box.
 */
      FixedVector<NDIM, T> getUpper() const { return upper; }

/**
 * Check if point is inside region.
 *
 * @param point Point to check.
 */
      bool isInside(T point[NDIM]) const;

/**
 * Extend the region by specified size on lower and upper sides of
 * region. This region is not modified but a new region is returned.
 *
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return extended region.
 */
      Region<NDIM, T> extend(T lowerExt[NDIM], T upperExt[NDIM]) const;

    private:
/** Lower and upper coordinates of region */
      Lucee::FixedVector<NDIM, T> lower, upper;
/** Volume of region */
      T volume;
  };

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(T shape[NDIM])
    : lower(0), upper(shape)
  {
    volume = 1;
    for (unsigned i=0; i<NDIM; ++i)
      volume *= shape[i];
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(T lo[NDIM], T up[NDIM])
    : lower(lo), upper(up)
  {
    volume = 1;
    for (unsigned i=0; i<NDIM; ++i)
      volume *= (up[i]-lo[i]);
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(const Region<NDIM, T>& rgn)
    : lower(rgn.lower), upper(rgn.upper), volume(rgn.volume)
  {
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>&
  Region<NDIM, T>::operator=(const Region<NDIM, T>& rgn)
  {
    if (&rgn == this)
      return *this;    
    lower = rgn.lower;
    upper = rgn.upper;
    volume = rgn.volume;
    return *this;
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

  template <unsigned NDIM, typename T>
  Region<NDIM, T>
  Region<NDIM, T>::extend(T lowerExt[NDIM], T upperExt[NDIM]) const
  {
  }
}

#endif // LC_REGION_H

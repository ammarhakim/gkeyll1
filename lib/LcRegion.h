/**
 * @file	LcRegion.h
 *
 * @brief	Region in an N-dimensional space.
 */

#ifndef LC_REGION_H
#define LC_REGION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
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
 * Creates an empty box, i.e. with 0 volume and same lower and upper
 * bounds.
 */
      Region();

/**
 * Create a region object with given shape. Lower limits are assumed
 * to be zero.
 *
 * @param shape Shape of region.
 */
      Region(const T shape[NDIM]);

/**
 * Create a region object with given lower and upper bounds.
 *
 * @param lower Lower limits of region.
 * @param upper Upper limits of region.
 */
      Region(const T lower[NDIM], const T upper[NDIM]);

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
 * Compare supplied region with this region.
 *
 * @param rgn Region object to compare.
 * @return true if regions are identical, false otherwise.
 */
      bool operator==(const Region<NDIM, T>& rgn) const;

/**
 * Compare supplied region with this region.
 *
 * @param rgn Region object to compare.
 * @return false if regions are identical, false otherwise.
 */
      bool operator!=(const Region<NDIM, T>& rgn) const;

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
 * @return true if point is inside region.
 */
      bool isInside(const T point[NDIM]) const;

/**
 * Return true if region is empty (zero or negative volume).
 *
 * @return True if region is empty.
 */
      bool isEmpty() const;

/**
 * Return intersection of supplied region with this one.
 *
 * @param rgn Region to intersect.
 * @return Intersected region.
 */
      Region<NDIM, T> intersect(const Region<NDIM, T>& rgn) const;

/**
 * Check if supplied region is inside this region. Tests true for self
 * containment.
 *
 * @param rgn Region to check.
 * @return True if this box contains rgn, false otherwise.
 */
      bool contains(const Region<NDIM, T>& rgn) const;

/**
 * Set lower bound of region in a given direction.
 *
 * @param dir Direction to reset bound.
 * @param lo New lower bound of region.
 */
      void setLower(unsigned dir, const T& lo) { 
        lower[dir] = lo;
        volume = calcVolume(); 
      }

/**
 * Set upper bound of region in a given direction.
 *
 * @param dir Direction to reset bound.
 * @param up New uppper bound of region.
 */
      void setUpper(unsigned dir, const T& up) {
        upper[dir] = up;
        volume = calcVolume(); 
      }

/**
 * Reset the bounds of box in specified direction to range
 * [nl,nu). This region is not modified but a new region is returned.
 *
 * @param dir Direction to reset bounds for.
 * @param nl New lower bound.
 * @param nu New upper bound.
 * @return region with new bounds
 */
      Region<NDIM, T> resetBounds(unsigned dir, const T nl, const T nu) const;

/**
 * Extend the region by specified size on lower and upper sides of
 * region. This region is not modified but a new region is returned.
 *
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return extended region.
 */
      Region<NDIM, T> extend(const T lowerExt[NDIM], const T upperExt[NDIM]) const;

/**
 * Inflate the region by adding one additional dimension by setting
 * lower and upper bounds of the dimension. The dimension is added as
 * the final dimension. A new region object is returned.
 *
 * @param lo Lower bound for new dimension.
 * @param up Upper bound for new dimension.
 * @return inflated region.
 */
      Region<NDIM+1, T> inflate(T lo, T up) const;

/**
 * Creates a new regions that is same as this one, except that
 * getUpper(dir) = getLower(dir)+1. In essence, the shape along 'dir'
 * is reduced to 1.
 *
 * @param dir Direction to deflate.
 * @return deflated region.
 */
      Region<NDIM, T> deflate(unsigned dir) const;

    private:
/** Lower coordinates of region */
      Lucee::FixedVector<NDIM, T> lower;
/** Upper coordinates of region */
      Lucee::FixedVector<NDIM, T> upper;
/** Volume of region */
      T volume;

/**
 * Compute volume of region.
 *
 * @return volume of region.
 */
      T calcVolume() const;
  };

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(const T shape[NDIM])
    : lower((T) 0), upper(shape)
  {
    volume = calcVolume();
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region(const T lo[NDIM], const T up[NDIM])
    : lower(lo), upper(up)
  {
    volume = calcVolume();
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
  Region<NDIM, T>::operator==(const Region<NDIM, T>& rgn) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      if ( (lower[i] != rgn.lower[i]) || (upper[i] != rgn.upper[i]) )
        return false;
    return true;
  }

  template <unsigned NDIM, typename T>
  bool
  Region<NDIM, T>::operator!=(const Region<NDIM, T>& rgn) const
  {
    return !this->operator==(rgn);
  }

  template <unsigned NDIM, typename T>
  bool
  Region<NDIM, T>::isInside(const T point[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      if ((point[i] < lower[i]) || (point[i] >= upper[i]))
        return false;
    return true;
  }

  template <unsigned NDIM, typename T>
  bool
  Region<NDIM, T>::isEmpty() const
  {
    for (unsigned i=0; i<NDIM; ++i)
      if (upper[i] <= lower[i])
        return true;
    return false;
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>
  Region<NDIM, T>::intersect(const Region<NDIM, T>& rgn) const
  {
    T lo[NDIM], up[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
// lower is max of regions' lower coordinates
      lo[i] = lower[i] > rgn.lower[i] ? lower[i] : rgn.lower[i];
// upper is min of regions' upper coordinates
      up[i] = upper[i] < rgn.upper[i] ? upper[i] : rgn.upper[i];
      if (up[i] <= lo[i])
        return Region<NDIM, T>();
    }
    return Region<NDIM, T>(lo, up);
  }

  template <unsigned NDIM, typename T>
  bool
  Region<NDIM, T>::contains(const Region<NDIM, T>& rgn) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      if ((rgn.lower[i] < lower[i]) || (rgn.lower[i] >= upper[i]))
        return false;
    for (unsigned i=0; i<NDIM; ++i)
      if ((rgn.upper[i] < lower[i]) || (rgn.upper[i] > upper[i]))
        return false;
    return true;
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>
  Region<NDIM, T>::resetBounds(unsigned dir, const T nl, const T nu) const
  {
    int newLo[NDIM], newUp[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      newLo[i] = lower[i];
      newUp[i] = upper[i];
    }
    newLo[dir] = nl;
    newUp[dir] = nu;
    return Region<NDIM, T>(newLo, newUp);
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>
  Region<NDIM, T>::extend(const T lowerExt[NDIM], const T upperExt[NDIM]) const
  {
    int newLo[NDIM], newUp[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      newLo[i] = lower[i]-lowerExt[i];
      newUp[i] = upper[i]+upperExt[i];
    }
    return Region<NDIM, T>(newLo, newUp);
  }

  template <unsigned NDIM, typename T>
  Region<NDIM+1, T>
  Region<NDIM, T>::inflate(T lo, T up) const
  {
    T newLo[NDIM+1], newUp[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i)
    {
      newLo[i] = lower[i];
      newUp[i] = upper[i];
    }
    newLo[NDIM] = lo;
    newUp[NDIM] = up;
    return Region<NDIM+1, T>(newLo, newUp);
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>
  Region<NDIM, T>::deflate(unsigned dir) const
  {
    if (dir >= NDIM)
    {
      Lucee::Except lce("Region::deflate: Deflate direction should be less than '");
      lce << NDIM << "'";
      throw lce;
    }
    return this->resetBounds(dir, lower[dir], lower[dir]+1);
  }

  template <unsigned NDIM, typename T>
  Region<NDIM, T>::Region()
    : lower(0), upper(0), volume(0)
  {
  }

  template <unsigned NDIM, typename T>
  T
  Region<NDIM, T>::calcVolume() const
  {
    T vol = 1;
    for (unsigned i=0; i<NDIM; ++i)
      vol *= getShape(i);
    return vol;
  }

/**
 * Create a new region from start and shape arrays.
 *
 * @param start Start indices of region.
 * @param shape Shape of region.
 * @return region with given start and shape.
 */
  template <unsigned NDIM, typename T>
  Region<NDIM, T>
  createRegionFromStartAndShape(const T start[NDIM], const T shape[NDIM])
  {
    T upper[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      upper[i] = start[i] + shape[i];
    return Region<NDIM, T>(start, upper);
  }
  
}

#endif // LC_REGION_H

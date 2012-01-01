/**
 * @file	LcDecompRegion.h
 *
 * @brief	A region that is divided into non-overlapping regions
 */

#ifndef LC_DECOMP_REGION_H
#define LC_DECOMP_REGION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcRegion.h>

// std includes
#include <vector>

namespace Lucee
{
  template <unsigned NDIM> 
  class DecompRegion
  {
    public:
/**
 * Create a new decompRegion object which decomposes the supplied
 * global region. The object created by this ctor creates a unitary
 * decomp: i.e. the global region is itself the sub-region. To perform
 * the actual decomposition use a derived class of
 * DecompRegionCalcIfc.
 *
 * @param globalRgn Global region to decompose.
 */
      DecompRegion(const Lucee::Region<NDIM, int>& globalRgn);

/**
 * Return bumber of sub-regions in decomposition.
 *
 * @return number of sub-regions.
 */
      unsigned getNumRegions() const;

/**
 * Return specified sub-region in decomposition.
 *
 * @param rn Region number.
 * @return region requested.
 */
      Lucee::Region<NDIM, int> getRegion(unsigned rn) const;

/**
 * Return specified sub-region in decomposition.
 *
 * @param rn Region number.
 * @return region requested.
 */
      Lucee::Region<NDIM, int> getGlobalRegion() const;

    private:
/** Global region */
      Lucee::Region<NDIM, int> globalRgn;
/** Regions making up decomposition */
      std::vector<Lucee::Region<NDIM, int> > rgns;
  };
}

#endif // LC_DECOMP_REGION_H

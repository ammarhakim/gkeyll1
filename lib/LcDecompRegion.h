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
/**
 * A decomposed region is an abstraction that holds a decomposition of
 * a given global region as (potentially) smaller non-overlapping
 * regions. The union of all the regions in the decomposition covers
 * the global region.
 */
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

/**
 * Get neigbors of target region for a given number of ghost cells on
 * each side of the target region. The neigbors also include corner
 * cells.
 *
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of neigbors region numbers.
 */
      std::vector<unsigned> getNeighbors(unsigned rn,
        const int lowerExt[NDIM], const int upperExt[NDIM]) const;

/**
 * Get neigbors of target region for a given number of ghost cells on
 * each side of the target region. The neigbors do not include corner
 * neigbors but only face neigbors.
 *
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of neigbors region numbers.
 */
      std::vector<unsigned> getFaceNeighbors(unsigned rn,
        const int lowerExt[NDIM], const int upperExt[NDIM]) const;

/**
 * Check if the decomposition covers global region. Returns true if it
 * does, false otherwise.
 *
 * @return return true if decomposition covers global region.
 */
      bool checkCovering() const;

    private:
/** Global region */
      Lucee::Region<NDIM, int> globalRgn;
/** Regions making up decomposition */
      std::vector<Lucee::Region<NDIM, int> > rgns;

/**
 * Clear decomposition to create new decomposition.
 */
      void clearDecomp();
  };
}

#endif // LC_DECOMP_REGION_H

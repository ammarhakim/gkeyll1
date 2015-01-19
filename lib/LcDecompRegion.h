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
#include <LcFixedVector.h>
#include <LcRegion.h>

// std includes
#include <functional>
#include <vector>
#include <map>

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
// forward declare algorithm calculation class so it can access our privates
      template <unsigned NNDIM> friend class DecompRegionCalcIfc;

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
 * Create a new decompRegion object from supplied one.
 *
 * @param decompRgn Decomposed region.
 */
      DecompRegion(const DecompRegion<NDIM>& decompRgn);

/**
 * Return number of sub-regions in decomposition.
 *
 * @return number of sub-regions.
 */
      unsigned getNumRegions() const;

/**
 * Return total number of sub-regions in decomposition, including
 * pseudo-regions.
 *
 * @return number of sub-regions.
 */
      unsigned getNumTotalRegions() const;

/**
 * Return specified sub-region in decomposition.
 *
 * @param rn Region number.
 * @return region requested.
 */
      Lucee::Region<NDIM, int> getRegion(unsigned rn) const;

/**
 * Return rank of specified sub-region in decomposition.
 *
 * @param rn Region number.
 * @return rank of region.
 */
      int getRank(unsigned rn) const;

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
 * cells. This method returns those neighbors from which we should
 * receive data.
 *
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of receive neigbors region numbers.
 */
      std::vector<unsigned> getRecvNeighbors(unsigned rn,
        const int lowerExt[NDIM], const int upperExt[NDIM]) const;

/**
 * Get neigbors of target region for a given number of ghost cells on
 * each side of the target region. The neigbors also include corner
 * cells. This method returns those neighbors to which we should
 * send data.
 *
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of send neigbors region numbers.
 */
      std::vector<unsigned> getSendNeighbors(unsigned rn,
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
 * Get list of regions that intersect the specified box.
 *
 * @param box Box to test intersection with.
 * @return List of regions that intersect box.
 */
      std::vector<unsigned> getIntersectingRegions(const Lucee::Region<NDIM, int> box) const;

/**
 * Check if the decomposition covers global region. Returns true if it
 * does, false otherwise.
 *
 * @return return true if decomposition covers global region.
 */
      bool checkCovering() const;

/**
 * Calculate ratio of minimum volume region to maximum volume region
 *
 * @return ratio of minimum volume to maximum volume region.
 */
      double calcMinMaxVolRatio() const;

/**
 * Compare if the supplied decomposed region is identical to this one.
 *
 * @param decompRgn Decomposed region to compare to.
 * @return true it regions are identical, false if not.
 */
      bool compareDecomp(const Lucee::DecompRegion<NDIM>& decompRgn) const;

    private:
/** Global region */
      Lucee::Region<NDIM, int> globalRgn;
/** Number of regions */
      unsigned numRegions;
/** Regions making up decomposition */
      std::vector<Lucee::Region<NDIM, int> > rgns;
/** Region to rank list */
      std::vector<int> rgnRank;

/**
 * A private class to define a comparison operator that allows using
 * ghost cell distributions as keys in a map. This is needed to cache
 * neighbor calculations that would otherwise be very expensive.
 */
      template <unsigned CDIM>
      struct FixedVecCmp : public std::binary_function<Lucee::FixedVector<CDIM, int>,
        Lucee::FixedVector<CDIM, int>, bool>
      {
/**
 * Function that compares two fixed vectors lexicographically.
 *
 * @param lhs Vector on left of < operator.
 * @param rhs Vector on right of < operator.
 * @return true is lhs < rhs, false otherwise
 */
          bool
          operator() (const Lucee::FixedVector<CDIM, int>& lhs,
            const Lucee::FixedVector<CDIM, int>& rhs) const
          {
            for (unsigned i=0; i<CDIM; ++i)
            {
              if (lhs[i] < rhs[i])
                return true;
              else if (lhs[i] > rhs[i])
                return false;
            }
            return false; // lhs == rhs hence lhs < rhs is false
          }
      };

/** Typedef map from ghost cell distributions -> neighbor list */
      typedef std::map<Lucee::FixedVector<2*NDIM, int>, std::vector<unsigned>, FixedVecCmp<2*NDIM> > NeighborMap_t;
/** Typedef pair of ghost cell distributions & neighbor list */
      typedef std::pair<Lucee::FixedVector<2*NDIM, int>, std::vector<unsigned> > NeighborPair_t;

/**
 * Private class to hold all neighbor information for each sub-region.
 */
      struct NeighborData
      {
/** Map to store ghost cell distribution -> neighbor list */
          mutable NeighborMap_t neighborMap;
      };

/** Map of sub-region number -> neighbor information for recv neighbors */
      mutable std::map<unsigned, NeighborData> rgnRecvNeighborMap;
/** Map of sub-region number -> neighbor information for send neighbors */
      mutable std::map<unsigned, NeighborData> rgnSendNeighborMap;

/**
 * Clear current decomposition to create new decomposition.
 */
      void clearDecomp();

/**
 * Append sub-region to decomposition. (This method is private but
 * DecompRegionCalcIfc class which is a friend can use it while
 * creating the decomposition).
 *
 * @param subRgn Region to add.
 */
      void addRegion(const Lucee::Region<NDIM, int>& subRgn);

/**
 * Append pseudo-region to decomposition. Pseudo regions are regions
 * created explicitly to handle periodic BCs. They are not used
 * otherwise.
 *
 * @param rank Rank of pseudo region.
 * @param subRgn Region to add.
 */
      void addPseudoRegion(int rank, const Lucee::Region<NDIM, int>& subRgn);

/**
 * Calculate neighbor information given ghost cell distribution. This
 * method returns those neighbors from which we should receive data.
 *
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of neigbors region numbers.
 */
      std::vector<unsigned> calcRecvNeighbors(unsigned rn,
        const int lowerExt[NDIM], const int upperExt[NDIM]) const;

/**
 * Calculate neighbor information given ghost cell distribution. This
 * method returns those neighbors to which we should send data.
 *
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of neigbors region numbers.
 */
      std::vector<unsigned> calcSendNeighbors(unsigned rn,
        const int lowerExt[NDIM], const int upperExt[NDIM]) const;

/**
 * Get neigbors of target region for a given number of ghost cells on
 * each side of the target region. The neigbors also include corner
 * cells.
 *
 * @param recvNeigh True to get receive neighbors, false to get send neighbors
 * @param rn Target region number
 * @param lowerExt Lenght of extension along lower end in each direction.
 * @param upperExt Lenght of extension along upper end in each direction.
 * @return list of neigbors region numbers.
 */
      std::vector<unsigned> getNeighbors(bool recvNeigh,
        unsigned rn, const int lowerExt[NDIM], const int upperExt[NDIM]) const;

/**
 * Compare if the supplied neighbor maps are identical.
 *
 * @param rgnNeighborMap1 First neighbor region map to compare.
 * @param rgnNeighborMap2 Second neighbor region map to compare.
 * @return true if identical, false if not.
 */
      bool cmpRgnNeighMap(const std::map<unsigned, NeighborData>& rgnNeighborMap1,
        const std::map<unsigned, NeighborData>& rgnNeighborMap2) const;
  };
}

#endif // LC_DECOMP_REGION_H

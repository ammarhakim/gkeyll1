/**
 * @file	LcMultiRegion.h
 *
 * @brief	A set of regions connected to each other.
 */

#ifndef LC_MULTI_REGION_H
#define LC_MULTI_REGION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee include
#include <LcRegion.h>

// std includes
#include <map>

namespace Lucee
{
/** Forward declare multi-region class so we can friend it */
  template <unsigned NDIM, typename T> class MultiRegion;

/** Enum to set lower/upper side */  
  enum RegionSide { LOWER, UPPER };

/**
 * This structure represents connectivity information for one edge of
 * a region. It is created by specifying the index of the target
 * region and the information (side and direction) of the connection
 * to the target region.
 */
  struct MultiRegionConnectivity
  {
    public:
/** Friend multi-region class so it can access our privates */
      template <unsigned NDIM, typename T> friend class MultiRegion;

/**
 * Create new default connectivity object for an edge. This ctor
 * creates an unconnected object.
 */
      MultiRegionConnectivity();

/**
 * Create new connectivity object for an edge of a region. If the edge
 * is not connected to any other region use -1 as the 'rgnIdx'. In
 * this case the other two parameters can be arbitrary.
 *
 * @param rgnIdx Index of target region. Use -1 to indicate edge is unconnected.
 * @param targetDir Direction of target region to which this edge is connected.
 * @param targetSide Use LOWER for lower side and UPPER for upper side.
 */
      MultiRegionConnectivity(int rgnIdx, unsigned targetDir, RegionSide targetSide);

    private:
/** Region index */
      int regionIndex;
/** Target direction */
      unsigned targetDir;
/** Target side */
      RegionSide targetSide;
  };

/**
 * A multi-region represents a set of N-dimensional rectangular
 * regions connected to each other to form a larger set (a
 * multi-region) in N-dimensional space. The connectivity between the
 * regions specifies the topology of the multi-region but not always
 * the geometry.
 */
  template <unsigned NDIM, typename T>
  class MultiRegion
  {
    private:
/**
 * Private structure to hold connectivity of a single region.
 */
      class RegionConn
      {
        public:
/**
 * Create new connectivity object. Connectivities are set to
 * unconnected state.
 *
 * @param rgn Region in multi-region object.
 */
          RegionConn(const Lucee::Region<NDIM, T>& rgn);

/** Region */
          Lucee::Region<NDIM, T> region;
/** Connectivity info of region connected to lower side */
          MultiRegionConnectivity lower[NDIM];
/** Connectivity info of region connected to upper side */
          MultiRegionConnectivity upper[NDIM];
      };

/**
 * Iterator to loop over regions and their neighbors.
 */
      class MultiRegionIterator
      {
        public:
/**
 * Create region iterator.
 *
 * @param mrPtr Pointer to multi-region object.
 */
          MultiRegionIterator(Lucee::MultiRegion<NDIM, T>& mrRef);

/** 
 * Increment iterator to next region. 
 */
          void next() { ++rgnMapItr; }

/**
 * Check if we are at the end.
 *
 * @return true if we are done, false otherwise.
 */
          bool atEnd() const { return rgnMapItr == multiRgnPtr->regionMap.end(); }

        private:
/** Pointer multi-region object */
          Lucee::MultiRegion<NDIM, T> *multiRgnPtr;
/** Region map iterator */
          typename std::map<int, RegionConn>::const_iterator rgnMapItr;
      };

    public:
/** Iterator for looping over region and neighbor data */
      typedef MultiRegionIterator iterator;

/**
 * Add a region with specified index. This region is assumed to be
 * unconnected. After all regions have been added, connectivities can
 * be specified using `setRegionLowerConnection` and
 * `setRegionUpperConnection` methods.
 *
 * @param idx Index of region.
 * @param rgn Region to add.
 */
      void addRegion(unsigned idx, const Lucee::Region<NDIM, T>& rgn);

/**
 * Number of regions in multi-region.
 *
 * @return number of regions in multi-region.
 */
      unsigned getNumRegions() const { return regionMap.size(); }

/**
 * Checks consistency of the multi-region. A multi-region is said to
 * be consistent if edges common to two regions have the same
 * shape. Note that a multi-region need not be consistent. In such
 * cases it is upto the user of this object to handle the
 * inconsistently connected regions.
 *
 * @return true if regions are consistent, false otherwise.
 */
      bool isConsistent() const;

 private:
/** Map of region ID -> regions */
      std::map<int, RegionConn> regionMap;
  };
}

#endif // LC_MULTI_REGION_H

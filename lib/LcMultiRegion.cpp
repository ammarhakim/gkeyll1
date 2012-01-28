/**
 * @file	LcMultiRegion.cpp
 *
 * @brief	A set of regions connected to each other.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee include
#include <LcMultiRegion.h>

// std includes
#include <vector>

namespace Lucee
{
  MultiRegionConnectivity::MultiRegionConnectivity()
    : regionIndex(-1), targetDir(0), targetSide(LOWER)
  {
  }

  MultiRegionConnectivity::MultiRegionConnectivity(
    int rgnIdx, unsigned targetDir, RegionSide targetSide)
    : regionIndex(rgnIdx), targetDir(targetDir), targetSide(targetSide)
  {
  }

  void
  MultiRegionConnectivity::reset(int rgnIdx, unsigned tD, RegionSide tS)
  {
    regionIndex = rgnIdx;
    targetDir = tD;
    targetSide = tS;
  }

  template <unsigned NDIM, typename T>
  MultiRegion<NDIM, T>::RegionConn::RegionConn(const Lucee::Region<NDIM, T>& rgn)
    : region(rgn)
  {
    for (unsigned i=0; i<NDIM; ++i)
    { // region is not connected
      lower[i] = MultiRegionConnectivity();
      upper[i] = MultiRegionConnectivity();
    }
  }

  template <unsigned NDIM, typename T>
  MultiRegion<NDIM, T>::MultiRegionIterator::MultiRegionIterator(Lucee::MultiRegion<NDIM, T>& mrRef)
    : multiRgnPtr(&mrRef)
  {
    rgnMapItr = multiRgnPtr->regionMap.begin();
  }

  template <unsigned NDIM, typename T>
  void
  MultiRegion<NDIM, T>::addRegion(unsigned idx, const Lucee::Region<NDIM, T>& rgn)
  {
// check if region with this ID already exists
    typename std::map<int, RegionConn>::iterator itr = regionMap.find(idx);
    if (itr == regionMap.end())
// add region information with unconnected neighbors
      regionMap.insert(std::pair<int, RegionConn>(idx, RegionConn(rgn)));
  }

// instantiations
  template class MultiRegion<1, float>;
  template class MultiRegion<2, float>;
  template class MultiRegion<3, float>;

  template class MultiRegion<1, double>;
  template class MultiRegion<2, double>;
  template class MultiRegion<3, double>;

  template class MultiRegion<1, int>;
  template class MultiRegion<2, int>;
  template class MultiRegion<3, int>;
}

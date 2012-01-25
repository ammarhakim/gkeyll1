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
    : regionIndex(-1), targetDir(0), targetSide(0)
  {
  }

  MultiRegionConnectivity::MultiRegionConnectivity(
    int rgnIdx, unsigned targetDir, unsigned targetSide)
    : regionIndex(rgnIdx), targetDir(targetDir), targetSide(targetSide)
  {
    if (rgnIdx != -1 && targetSide > 1)
    {
      Lucee::Except lce("MultiRegionConnectivity::MultiRegionConnectivity: ");
      lce << " For connected edge side can be 0 or 1. Specified " << targetSide
          << " instead";
      throw lce;
    }
  }

  void
  MultiRegionConnectivity::reset(int rgnIdx, unsigned tD, unsigned tS)
  {
    if (rgnIdx != -1 && tS > 1)
    {
      Lucee::Except lce("MultiRegionConnectivity::reset: ");
      lce << " For connected edge side can be 0 or 1. Specified " << tS
          << " instead";
      throw lce;
    }

    regionIndex = rgnIdx;
    targetDir = tD;
    targetSide = tS;
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

/**
 * @file	LcDecompRegion.cpp
 *
 * @brief	A region that is divided into non-overlapping regions
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcArray.h>
#include <LcDecompRegion.h>

// std includes
#include <limits>

namespace Lucee
{
  template <unsigned NDIM> 
  DecompRegion<NDIM>::DecompRegion(const Lucee::Region<NDIM, int>& globalRgn)
    : globalRgn(globalRgn)
  {
    rgns.push_back(globalRgn); // by default global region is not decomposed
  }

  template <unsigned NDIM> 
  unsigned
  DecompRegion<NDIM>::getNumRegions() const
  {
    return rgns.size();
  }

  template <unsigned NDIM> 
  Lucee::Region<NDIM, int>
  DecompRegion<NDIM>::getRegion(unsigned rn) const
  {
    if (rn >= rgns.size())
    {
      Lucee::Except lce("DecompRegion::getRegion: Region number ");
      lce << rn << " is not valid. Should be smaller than " << rgns.size();
      throw lce;
    }
    return rgns[rn];
  }

  template <unsigned NDIM> 
  Lucee::Region<NDIM, int>
  DecompRegion<NDIM>::getGlobalRegion() const
  {
    return globalRgn;
  }

  template <unsigned NDIM> 
  std::vector<unsigned>
  DecompRegion<NDIM>::getNeighbors(unsigned rn,
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
// create vector to identify ghost cell distribution
    Lucee::FixedVector<2*NDIM, int> gcd(0);
    for (unsigned i=0; i<NDIM; ++i)
    {
      gcd[i] = lowerExt[i];
      gcd[NDIM+i] = upperExt[i];
    }

// first check if any data exists for this region number
    typename std::map<unsigned, NeighborData>::const_iterator rgnItr
      = rgnNeighborMap.find(rn);
    if (rgnItr != rgnNeighborMap.end())
    {
// check if this ghost cell distribution is computed
      typename NeighborMap_t::const_iterator gstItr
        = rgnItr->second.neighborMap.find(gcd);
      if (gstItr != rgnItr->second.neighborMap.end())
        return gstItr->second;
    }
// at this point we need to compute neighbors
    std::vector<unsigned> ninfo = calcNeighbors(rn, lowerExt, upperExt);
// ADD THIS
  }
  
  template <unsigned NDIM> 
  void
  DecompRegion<NDIM>::clearDecomp()
  {
    rgns.clear();
  }

  template <unsigned NDIM> 
  bool
  DecompRegion<NDIM>::checkCovering() const
  {
// create array over global region
    Lucee::Array<NDIM, int> check(globalRgn, 0);
// loop over each sub-region
    typename std::vector<Lucee::Region<NDIM, int> >::const_iterator itr
      = rgns.begin();
    for ( ; itr != rgns.end(); ++itr)
    {
      Lucee::ColMajorSequencer<NDIM> seq(*itr);
// loop over region, incrementing count in 'check' array
      while (seq.step())
        check(seq.getIndex()) += 1;
    }
// if any location is visited more than once (or not at all) it will
// have a number other than 1 in it.
    Lucee::ColMajorSequencer<NDIM> seq(globalRgn);
    while (seq.step())
    {
      if (check(seq.getIndex()) != 1)
        return false;
    }
    return true;
  }

  template <unsigned NDIM> 
  double
  DecompRegion<NDIM>::calcMinMaxVolRatio() const
  {
    int minVol = std::numeric_limits<int>::max();
    int maxVol = 0;
// loop over each sub-region
    typename std::vector<Lucee::Region<NDIM, int> >::const_iterator itr
      = rgns.begin();
    for ( ; itr != rgns.end(); ++itr)
    {
      int vol = itr->getVolume();
      minVol = vol < minVol ? vol : minVol;
      maxVol = vol > maxVol ? vol : maxVol;
    }
    return (double) minVol / (double) maxVol;
  }

  template <unsigned NDIM> 
  void
  DecompRegion<NDIM>::addRegion(const Lucee::Region<NDIM, int>& subRgn)
  {
    rgns.push_back(subRgn);
  }

  template <unsigned NDIM> 
  std::vector<unsigned> 
  DecompRegion<NDIM>::calcNeighbors(unsigned rn,
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    std::vector<unsigned> nl;
// get current region 
    Lucee::Region<NDIM, int> currRgn = getRegion(rn);
// intersect it with all other regions
    for (unsigned i=0; i<getNumRegions(); ++i)
    {
      if (rn == i) 
        continue; // no need to intersect with ourself
// check intersection
      if ( ! currRgn.intersect( getRegion(i) ).isEmpty() )
        nl.push_back(i);
    }
    return nl;
  }

// instantiations
  template class DecompRegion<1>;
  template class DecompRegion<2>;
  template class DecompRegion<3>;
  template class DecompRegion<4>;
  template class DecompRegion<5>;
  template class DecompRegion<6>;
  template class DecompRegion<7>;
}

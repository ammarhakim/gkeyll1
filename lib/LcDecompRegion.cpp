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
  DecompRegion<NDIM>::DecompRegion(const DecompRegion<NDIM>& decompRgn)
    : globalRgn(decompRgn.globalRgn), rgns(decompRgn.rgns),
      rgnNeighborMap(decompRgn.rgnNeighborMap)
  {
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
// NOTE: This method checks if the neighbor calculation for this
// sub-region and ghost cell distribution was requested before. If so
// it just returns the previously computed value. Otherwise, it
// computes the neighbors, stores them in the appropriate map and
// returns the newly computed neighbors. The neighbor calculation is
// quite expensive, specially in 2D and 3D so the caching can improve
// performance. The code below is complicated by the fact that for
// each sub-region, neighbor information for different ghost
// distribution might be requested. This requires the use of a
// "two-layer" map.

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
    { // region exists
// check if this ghost cell distribution is computed
      typename NeighborMap_t::const_iterator gstItr
        = rgnItr->second.neighborMap.find(gcd);
      if (gstItr != rgnItr->second.neighborMap.end())
      { // ghost cell distribution exists
        return gstItr->second;
      }
      else
      { // ghost cell distribution does not exist
// compute neighbors
        std::vector<unsigned> ninfo = calcNeighbors(rn, lowerExt, upperExt);
// insert into map so next time we need not compute neighbors all over again
        rgnItr->second.neighborMap.insert(NeighborPair_t(gcd, ninfo));
        return ninfo;
      }
    }
    else
    { // region does not exist
// compute neighbors
        std::vector<unsigned> ninfo = calcNeighbors(rn, lowerExt, upperExt);
// create new neighbor data object
        NeighborData ndat;
        ndat.neighborMap.insert(NeighborPair_t(gcd, ninfo));
// now insert this into region map
        rgnNeighborMap.insert(
          std::pair<unsigned, NeighborData>(rn, ndat));
        return ninfo;
    }
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
  bool
  DecompRegion<NDIM>::compareDecomp(const Lucee::DecompRegion<NDIM>& decompRgn) const
  {
    if (globalRgn != decompRgn.globalRgn)
      return false;

    if (rgns.size() != decompRgn.rgns.size()) 
      return false;
    for (unsigned i=0; i<rgns.size(); ++i)
      if (rgns[i] != decompRgn.rgns[i])
        return false;

    if (rgnNeighborMap.size() != decompRgn.rgnNeighborMap.size())
      return false;

    typename std::map<unsigned, NeighborData>::const_iterator rnmItr
      = rgnNeighborMap.begin();
    for ( ; rnmItr != rgnNeighborMap.end(); ++rnmItr)
    {
      typename std::map<unsigned, NeighborData>::const_iterator d_rnmItr
        = decompRgn.rgnNeighborMap.find(rnmItr->first);
      if (d_rnmItr == decompRgn.rgnNeighborMap.end()) 
        return false;

      if (rnmItr->second.neighborMap.size() != d_rnmItr->second.neighborMap.size())
        return false;

      typename NeighborMap_t::const_iterator nmItr = 
        rnmItr->second.neighborMap.begin();
      for ( ; nmItr != rnmItr->second.neighborMap.end(); ++nmItr)
      {
        typename NeighborMap_t::const_iterator d_nmItr
          = d_rnmItr->second.neighborMap.find(nmItr->first);
        if (d_nmItr == d_rnmItr->second.neighborMap.end()) 
          return false;

        if (nmItr->second.size() != d_nmItr->second.size()) 
          return false;

        for (unsigned k=0; k<nmItr->second.size(); ++k)
          if (nmItr->second[k] != d_nmItr->second[k])
            return false;
      }
    }

    return true;
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
      if ( ! currRgn.extend(lowerExt, upperExt).intersect( getRegion(i) ).isEmpty() )
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

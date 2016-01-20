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
    : globalRgn(globalRgn), numRegions(0)
  {
    addRegion(globalRgn);
  }

  template <unsigned NDIM> 
  DecompRegion<NDIM>::DecompRegion(const DecompRegion<NDIM>& decompRgn)
    : globalRgn(decompRgn.globalRgn), numRegions(decompRgn.numRegions),
      rgns(decompRgn.rgns), rgnRank(decompRgn.rgnRank),
      rgnRecvNeighborMap(decompRgn.rgnRecvNeighborMap),
      rgnSendNeighborMap(decompRgn.rgnSendNeighborMap)
  {
  }

  template <unsigned NDIM> 
  unsigned
  DecompRegion<NDIM>::getNumRegions() const
  {
    return numRegions;
  }

  template <unsigned NDIM> 
  unsigned
  DecompRegion<NDIM>::getNumTotalRegions() const
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
  int
  DecompRegion<NDIM>::getRank(unsigned rn) const
  {
    if (rn >= rgnRank.size())
    {
      Lucee::Except lce("DecompRegion::getRank: Region number ");
      lce << rn << " is not valid. Should be smaller than " << rgns.size();
      throw lce;
    }
    return rgnRank[rn];
  }

  template <unsigned NDIM> 
  Lucee::Region<NDIM, int>
  DecompRegion<NDIM>::getGlobalRegion() const
  {
    return globalRgn;
  }

  template <unsigned NDIM>
  std::vector<unsigned>
  DecompRegion<NDIM>::getRecvNeighbors(unsigned rn,
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    return getNeighbors(true, rn, lowerExt, upperExt);
  }

  template <unsigned NDIM>
  std::vector<unsigned>
  DecompRegion<NDIM>::getSendNeighbors(unsigned rn,
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    return getNeighbors(false, rn, lowerExt, upperExt);
  }

  
  template <unsigned NDIM> 
  void
  DecompRegion<NDIM>::clearDecomp()
  {
    numRegions = 0;
    rgns.clear();
    rgnRank.clear();
  }

  template <unsigned NDIM> 
  std::vector<unsigned>
  DecompRegion<NDIM>::getIntersectingRegions(const Lucee::Region<NDIM, int> box) const
  {
    std::vector<unsigned> iRgns;
    for (unsigned i=0; i<getNumRegions(); ++i)
      if (box.intersect(rgns[i]).getVolume() > 0)
        iRgns.push_back(i);
    return iRgns;
  }

  template <unsigned NDIM> 
  bool
  DecompRegion<NDIM>::checkCovering() const
  {
// NEED TO FIX THIS AS THE ARRAY IS A MONSTER: ONE WAY AROUND THIS TO
// CONVERT THE CHECK TO PER-DOMAIN AND NOT PER-CELL

//    create array over global region
    Lucee::Array<NDIM, int> check(globalRgn, 0);
// loop over each sub-region
    for (unsigned i=0; i<getNumRegions(); ++i)
    {
      Lucee::ColMajorSequencer<NDIM> seq(rgns[i]);
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
    for (unsigned i=0; i<getNumRegions(); ++i)
    {
      int vol = rgns[i].getVolume();
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

    return cmpRgnNeighMap(rgnRecvNeighborMap, decompRgn.rgnRecvNeighborMap)
      && cmpRgnNeighMap(rgnSendNeighborMap, decompRgn.rgnSendNeighborMap);
  }

  template <unsigned NDIM> 
  void
  DecompRegion<NDIM>::addRegion(const Lucee::Region<NDIM, int>& subRgn)
  {
    rgns.push_back(subRgn);
    rgnRank.push_back(numRegions);
    numRegions++;
  }

  template <unsigned NDIM> 
  void
  DecompRegion<NDIM>::addPseudoRegion(int rank, const Lucee::Region<NDIM, int>& subRgn)
  {
    rgns.push_back(subRgn);
// Pseudo-regions will have ranks that do not correspond to their
// region number. So we need the rank specified explicitly.
    rgnRank.push_back(rank);
  }

  template <unsigned NDIM> 
  std::vector<unsigned> 
  DecompRegion<NDIM>::calcRecvNeighbors(unsigned rn,
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    std::vector<unsigned> nl;
// get current region 
    Lucee::Region<NDIM, int> currRgn = getRegion(rn);
// intersect extended region it with all other regions
    for (unsigned i=0; i<getNumTotalRegions(); ++i)
    {
      if (rn == i) 
        continue; // no need to intersect with ourself
// check intersection
      if ( ! currRgn.extend(lowerExt, upperExt).intersect( getRegion(i) ).isEmpty() )
        nl.push_back(i);
    }
    return nl;
  }

  template <unsigned NDIM> 
  std::vector<unsigned> 
  DecompRegion<NDIM>::calcSendNeighbors(unsigned rn,
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    std::vector<unsigned> nl;
// get current region 
    Lucee::Region<NDIM, int> currRgn = getRegion(rn);
// intersect it with all other regions extended by ghost cells
    for (unsigned i=0; i<getNumTotalRegions(); ++i)
    {
      if (rn == i) 
        continue; // no need to intersect with ourself
// check intersection
      if ( ! currRgn.intersect( getRegion(i).extend(lowerExt, upperExt) ).isEmpty() )
        nl.push_back(i);
    }
    return nl;
  }

  template <unsigned NDIM> 
  std::vector<unsigned>
  DecompRegion<NDIM>::getNeighbors(bool recvNeigh, unsigned rn,
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

// set reference to correct map to use
    std::map<unsigned, NeighborData>& rnm = 
      recvNeigh ? rgnRecvNeighborMap : rgnSendNeighborMap;

// create vector to identify ghost cell distribution
    Lucee::FixedVector<2*NDIM, int> gcd(0);
    for (unsigned i=0; i<NDIM; ++i)
    {
      gcd[i] = lowerExt[i];
      gcd[NDIM+i] = upperExt[i];
    }

// first check if any data exists for this region number
    typename std::map<unsigned, NeighborData>::const_iterator rgnItr
      = rnm.find(rn);
    if (rgnItr != rnm.end())
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
        std::vector<unsigned> ninfo = recvNeigh ?
          calcRecvNeighbors(rn, lowerExt, upperExt) : calcSendNeighbors(rn, lowerExt, upperExt);
// insert into map so next time we need not compute neighbors all over again
        rgnItr->second.neighborMap.insert(NeighborPair_t(gcd, ninfo));
        return ninfo;
      }
    }
    else
    { // region does not exist
// compute neighbors
        std::vector<unsigned> ninfo = recvNeigh ?
          calcRecvNeighbors(rn, lowerExt, upperExt) : calcSendNeighbors(rn, lowerExt, upperExt);
// create new neighbor data object
        NeighborData ndat;
        ndat.neighborMap.insert(NeighborPair_t(gcd, ninfo));
// now insert this into region map
        rnm.insert(
          std::pair<unsigned, NeighborData>(rn, ndat));
        return ninfo;
    }
  }

  template <unsigned NDIM>
  bool
  DecompRegion<NDIM>::cmpRgnNeighMap(const std::map<unsigned, NeighborData>& rgnNeighborMap1,
    const std::map<unsigned, NeighborData>& rgnNeighborMap2) const
  {
    if (rgnNeighborMap1.size() != rgnNeighborMap2.size())
      return false;

    typename std::map<unsigned, NeighborData>::const_iterator rnmItr
      = rgnNeighborMap1.begin();
    for ( ; rnmItr != rgnNeighborMap1.end(); ++rnmItr)
    {
      typename std::map<unsigned, NeighborData>::const_iterator d_rnmItr
        = rgnNeighborMap2.find(rnmItr->first);
      if (d_rnmItr == rgnNeighborMap2.end()) 
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

// instantiations
  template class DecompRegion<1>;
  template class DecompRegion<2>;
  template class DecompRegion<3>;
  template class DecompRegion<4>;
  template class DecompRegion<5>;
  template class DecompRegion<6>;
  template class DecompRegion<7>;
}

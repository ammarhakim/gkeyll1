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
#include <LcDecompRegion.h>

namespace Lucee
{
  template <unsigned NDIM> 
  DecompRegion<NDIM>::DecompRegion(const Lucee::Region<NDIM, int>& globalRgn)
    : globalRgn(globalRgn)
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

// instantiations
  template class DecompRegion<1>;
  template class DecompRegion<2>;
  template class DecompRegion<3>;
  template class DecompRegion<4>;
  template class DecompRegion<5>;
  template class DecompRegion<6>;
  template class DecompRegion<7>;
}

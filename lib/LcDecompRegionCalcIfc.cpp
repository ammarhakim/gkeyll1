/**
 * @file	LcDecompRegionCalcIfc.cpp
 *
 * @brief	Base class for decomposition algorithms.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegionCalcIfc.h>
#include <LcRowMajorSequencer.h>

// std includes
#include <iostream>

namespace Lucee
{
// set module name
  template <> const char *DecompRegionCalcIfc<1>::id = "DecompRegionCalc1D";
  template <> const char *DecompRegionCalcIfc<2>::id = "DecompRegionCalc2D";
  template <> const char *DecompRegionCalcIfc<3>::id = "DecompRegionCalc3D";
  template <> const char *DecompRegionCalcIfc<4>::id = "DecompRegionCalc4D";
  template <> const char *DecompRegionCalcIfc<5>::id = "DecompRegionCalc5D";

  template <unsigned NDIM> 
  DecompRegionCalcIfc<NDIM>::DecompRegionCalcIfc()
  {
    for (unsigned d=0; d<NDIM; ++d)
      isPeriodic[d] = false;
  }

  template <unsigned NDIM> 
  void
  DecompRegionCalcIfc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
  }

  template <unsigned NDIM> 
  void
  DecompRegionCalcIfc<NDIM>::setPeriodicDir(unsigned dir, bool isp)
  {
    if (dir>=0 && dir<NDIM)
      isPeriodic[dir] = isp;
    else
    {
      Lucee::Except lce(
        "DecompRegionCalcIfc::setPeriodicDirs: Direction must be between 0 and ");
      lce << (NDIM-1);
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::calcDecomp(unsigned nrgns, Lucee::DecompRegion<NDIM>& decompRgn)
  {
    decompRgnPtr = &decompRgn; // so addRegion can access
// clear existing decomposition
    decompRgn.clearDecomp();
// call derived class to compute decomposition
    this->decompose(nrgns, decompRgn.getGlobalRegion());
// check if decompostion covers global region
    if (decompRgn.checkCovering() == false)
      throw Lucee::Except(
        "DecompRegionCalcIfc: Decomposition does not cover global region");

// take into account periodic directions (if any)
    handlePeriodicDirs();
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::addRegion(const Lucee::Region<NDIM, int>& subRgn)
  {
    decompRgnPtr->addRegion(subRgn);
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::handlePeriodicDirs()
  {
// The basic idea here is to duplicate the decomposition in the
// periodic directions. This is a bit tricky as periodicity of corners
// also needs to be constructed. So, for example, in 2D if the
// decomposition is unitary (i.e. only a single sub-region) and the
// domain is periodic in both directions, then EIGHT additional
// regions will be added. It is possible that for a large problem
// (1000s of regions) there is a huge amount of storage wasted for
// periodic BCs, but I think that is unavoidable, in general. However,
// it is possible that one could duplicate only the regions that lie
// on a boundary, but that might be a much more trickier
// implementation and will not work in general. (A. Hakim 2/06/2013)

// This algorithm was originally written by Mahmood Miah and A. Hakim
// for the Facets project. This implementation is essentially the same
// as in Facets. I am retaining the original comments. (A. Hakim,
// 2/5/2013)

    Lucee::Region<NDIM, int> rgn = decompRgnPtr->getGlobalRegion();
// make a copy as list of regions changes as it is extended.
    std::vector<Lucee::Region<NDIM, int> > realVec;
    for (unsigned i=0; i<decompRgnPtr->getNumRegions(); ++i)
      realVec.push_back( decompRgnPtr->getRegion(i) );

    int itrMin[NDIM], itrMax[NDIM];
    for (size_t i = 0; i < NDIM; ++i) {
      itrMin[i] = isPeriodic[i] ? -1 : 0;
      itrMax[i] = isPeriodic[i] ? 2 : 1;
    }
    Lucee::Region<NDIM, int> shift(itrMin, itrMax);
    Lucee::RowMajorSequencer<NDIM> seq(shift);
    size_t nTotBoxes = decompRgnPtr->getNumRegions();
    int idx[NDIM];
    while ( seq.step() ) {
      seq.fillWithIndex(idx);
// skip case of zero offset
      bool skipStep = true;
      for (size_t i=0; i<NDIM; ++i) {
        if (idx[i] != 0) {
          skipStep = false;
          break;
        }
      }
// shift boxes, renumber, and add to the box map
      if (!skipStep) {
        int lext[NDIM], uext[NDIM], dist[NDIM];
        for (size_t i = 0; i < NDIM; ++i) {
          dist[i] = rgn.getShape(i);
          lext[i] = -1*idx[i]*dist[i];
          uext[i] =  1*idx[i]*dist[i];
        }
        for (unsigned i=0 ; i<realVec.size(); ++i)
          decompRgnPtr->addPseudoRegion(i, realVec[i].extend(lext, uext));
      }
    }
  }

// instantiations
  template class DecompRegionCalcIfc<1>;
  template class DecompRegionCalcIfc<2>;
  template class DecompRegionCalcIfc<3>;
  template class DecompRegionCalcIfc<4>;
  template class DecompRegionCalcIfc<5>;
}

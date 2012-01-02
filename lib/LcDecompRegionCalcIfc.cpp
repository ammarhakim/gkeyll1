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

namespace Lucee
{
  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::calcDecomp(unsigned nrgns, Lucee::DecompRegion<NDIM>& decompRgn)
  {
    decompRgnPtr = &decompRgn; // set so addRegion can access
// clear existing decomposition
    decompRgn.clearDecomp();
// call derived class to compute decomposition
    this->decompose(nrgns, decompRgn.getGlobalRegion());
// check if decompostion covers global region
    if (decompRgn.checkCovering() == false)
      throw Lucee::Except(
        "DecompRegionCalcIfc: Decomposition does not cover global region");
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::addRegion(const Lucee::Region<NDIM, int>& subRgn)
  {
    decompRgnPtr->addRegion(subRgn);
  }

// instantiations
  template class DecompRegionCalcIfc<1>;
  template class DecompRegionCalcIfc<2>;
  template class DecompRegionCalcIfc<3>;
  template class DecompRegionCalcIfc<4>;
  template class DecompRegionCalcIfc<5>;
  template class DecompRegionCalcIfc<6>;
  template class DecompRegionCalcIfc<7>;
}

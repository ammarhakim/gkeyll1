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
// set module name
  template <> const char *DecompRegionCalcIfc<1>::id = "DecompRegionCalc1D";
  template <> const char *DecompRegionCalcIfc<2>::id = "DecompRegionCalc2D";
  template <> const char *DecompRegionCalcIfc<3>::id = "DecompRegionCalc3D";
  template <> const char *DecompRegionCalcIfc<4>::id = "DecompRegionCalc4D";
  template <> const char *DecompRegionCalcIfc<5>::id = "DecompRegionCalc5D";
  template <> const char *DecompRegionCalcIfc<6>::id = "DecompRegionCalc6D";
  template <> const char *DecompRegionCalcIfc<7>::id = "DecompRegionCalc7D";

  template <unsigned NDIM> 
  void
  DecompRegionCalcIfc<NDIM>::readInput(Lucee::LuaTable& tbl)
  { // does nothing
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

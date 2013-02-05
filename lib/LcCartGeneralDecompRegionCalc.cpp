/**
 * @file	LcCartGeneralDecompRegionCalc.cpp
 *
 * @brief	Decomposition with specified cuts.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartGeneralDecompRegionCalc.h>
#include <LcMatrix.h>

// std includes
#include <cmath>

namespace Lucee
{
// set constructor name
  template <> const char *CartGeneralDecompRegionCalc<1>::id = "CartGeneral";
  template <> const char *CartGeneralDecompRegionCalc<2>::id = "CartGeneral";
  template <> const char *CartGeneralDecompRegionCalc<3>::id = "CartGeneral";
  template <> const char *CartGeneralDecompRegionCalc<4>::id = "CartGeneral";

  template <unsigned NDIM>
  CartGeneralDecompRegionCalc<NDIM>::CartGeneralDecompRegionCalc()
    : nsub(1)
  {
  }

  template <unsigned NDIM>
  CartGeneralDecompRegionCalc<NDIM>::CartGeneralDecompRegionCalc(unsigned nrgns)
    : nsub(nrgns)
  {
  }

  template <unsigned NDIM>
  void
  CartGeneralDecompRegionCalc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    DecompRegionCalcIfc<NDIM>::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  CartGeneralDecompRegionCalc<NDIM>::decompose(unsigned nrgs, const Lucee::Region<NDIM, int>& globalRgn)
  {
// check if cuts match up with number of regions
    if (nrgs != nsub)
    {
      Lucee::Except lce("CartGeneralDecompRegionCalc::decompose: Number of sub-regions ");
      lce << nsub << ") does not match requested sub-regions (" << nrgs << ")";
      throw lce;
    }
  }

// instantiations
  template class CartGeneralDecompRegionCalc<1>;
  template class CartGeneralDecompRegionCalc<2>;
  template class CartGeneralDecompRegionCalc<3>;
  template class CartGeneralDecompRegionCalc<4>;
}

/**
 * @file	LcCartProdDecompRegionCalc.cpp
 *
 * @brief	Decomposition with specified cuts.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>

namespace Lucee
{
  template <unsigned NDIM>
  CartProdDecompRegionCalc<NDIM>::CartProdDecompRegionCalc(const unsigned c[NDIM])
  {
    nsub = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      cuts[i] = c[i];
      nsub *= cuts[i];
    }
  }

  template <unsigned NDIM>
  void
  CartProdDecompRegionCalc<NDIM>::decompose(unsigned nrgs, const Lucee::Region<NDIM, int>& globalRgn)
  {
// check if cuts match up with number of regions
    if (nrgs != nsub)
    {
      Lucee::Except lce("CartProdDecompRegionCalc::decompose: Number of sub-regions from cuts (");
      lce << nsub << ") does not match requested sub-regions (" << nrgs << ")";
      throw lce;
    }
  }

// instantiations
  template class CartProdDecompRegionCalc<1>;
  template class CartProdDecompRegionCalc<2>;
  template class CartProdDecompRegionCalc<3>;
  template class CartProdDecompRegionCalc<4>;
  template class CartProdDecompRegionCalc<5>;
  template class CartProdDecompRegionCalc<6>;
  template class CartProdDecompRegionCalc<7>;
}

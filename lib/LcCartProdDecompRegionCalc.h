/**
 * @file	LcCartProdDecompRegionCalc.h
 *
 * @brief	Decomposition with specified cuts.
 */

#ifndef LC_CART_PROD_DECOMP_REGION_CALC_H
#define LC_CART_PROD_DECOMP_REGION_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegionCalcIfc.h>

namespace Lucee
{
/**
 * Simple decomposition algorithm in which number of divisions in each
 * direction is explicitly specified.
 */
  template <unsigned NDIM>
  class CartProdDecompRegionCalc : public Lucee::DecompRegionCalcIfc<NDIM>
  {
    public:
/**
 * Create a new decomp calculator with specified number of cuts in
 * each direction. The value cut[n] is the number of divisions in the
 * n-th direction.
 *
 * @param cuts Value cuts[n] is number of divisions along n-th direction.
 */
      CartProdDecompRegionCalc(const unsigned cuts[NDIM]);

    protected:

/**
 * Method used to compute actual decomposition.
 *
 * @param nrgns Number of regions to create.
 * @param globalRgn Global region to decompose.
 */
      void decompose(unsigned nrgs, const Lucee::Region<NDIM, int>& globalRgn);

    private:
/** Number of cuts in each direction */
      unsigned cuts[NDIM];
/** Number of total sub-regions */
      unsigned nsub;
  };
}

#endif // LC_DECOMP_REGION_CALC_IFC_H

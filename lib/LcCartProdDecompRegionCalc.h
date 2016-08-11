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
#include <LcColMajorIndexer.h>
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
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Default constructor: needed to allow creation from Lua.
 */
      CartProdDecompRegionCalc();

/**
 * Create a new decomp calculator with specified number of cuts in
 * each direction. The value cut[n] is the number of divisions in the
 * n-th direction.
 *
 * @param cuts Value cuts[n] is number of divisions along n-th direction.
 */
      CartProdDecompRegionCalc(const unsigned cuts[NDIM]);

/**
 * Create a new decomp calculator with specified number of cuts in
 * each direction. The value cut[n] is the number of divisions in the
 * n-th direction.
 *
 * @param cuts Value cuts[n] is number of divisions along n-th direction.
 */
      CartProdDecompRegionCalc(const int cuts[NDIM]);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Get the cuts for this decomposition.
 *
 * @param cuts  On output, cuts used in this decomposition.
 */
      void fillWithCuts(int cuts[NDIM]) const;

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
      int cuts[NDIM];
/** Number of total sub-regions */
      unsigned nsub;
/** Indexer for mapping cut-index index to rank */
      Lucee::ColMajorIndexer<NDIM> cutIndexer;

/**
 * Set number cuts for use in decomposition.
 *
 * @param cuts Value cuts[n] is number of divisions along n-th direction.
 */
      void setCuts(const unsigned cuts[NDIM]);
  };
}

#endif // LC_CART_PROD_DECOMP_REGION_CALC_H

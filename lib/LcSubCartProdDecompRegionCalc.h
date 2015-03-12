/**
 * @file	LcSubCartProdDecompRegionCalc.h
 *
 * @brief	Decomposition based on higher-dimensional decomposition
 */

#ifndef LC_SUB_CART_PROD_DECOMP_REGION_CALC_H
#define LC_SUB_CART_PROD_DECOMP_REGION_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegionCalcIfc.h>

namespace Lucee
{
/**
 * This decomposition object is used a higher-dimensional
 * decomposition object to create a lower dimensional one. The HDIM
 * parameter is the dimension of the higher-dimension decomposition
 * object which should be used.
 */
  template <unsigned NDIM, unsigned HDIM>
  class SubCartProdDecompRegionCalc : public Lucee::DecompRegionCalcIfc<NDIM>
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Default constructor: needed to allow creation from Lua.
 */
      SubCartProdDecompRegionCalc();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

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

/**
 * Set number cuts for use in decomposition.
 *
 * @param cuts Value cuts[n] is number of divisions along n-th direction.
 */
      void setCuts(const int cuts[NDIM]);

/**
 * Set communicator valid for local rank from list of communicators.
 *
 * @param clist List of communicators
 */
      void setValidComm(const std::vector<TxCommBase*>& clist);
  };
}

#endif // LC_SUB_CART_PROD_DECOMP_REGION_CALC_H

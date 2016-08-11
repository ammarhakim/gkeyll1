/**
 * @file	LcCartGeneralDecompRegionCalc.h
 *
 * @brief	Decomposition into arbitrary number of regions.
 */

#ifndef LC_CART_GENERAL_DECOMP_REGION_CALC_H
#define LC_CART_GENERAL_DECOMP_REGION_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegionCalcIfc.h>

namespace Lucee
{
/**
 * This class decomposes a given region into an arbitrary number of
 * sub-regions.
 */
  template <unsigned NDIM>
  class CartGeneralDecompRegionCalc : public Lucee::DecompRegionCalcIfc<NDIM>
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Default constructor: needed to allow creation from Lua.
 */
      CartGeneralDecompRegionCalc();

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
/**
 * Creates a set of boxes
 *
 * @param breakme box which should be broken down
 * @param num number of boxes to put in each direction
 * @return a vector of decomposed FcBoxes
 */
      std::vector< Lucee::Region<NDIM, int> > breakBoxes(size_t subBox[NDIM], 
        const Lucee::Region<NDIM, int>& box);
  };
}

#endif // LC_CART_GENERAL_DECOMP_REGION_CALC_H

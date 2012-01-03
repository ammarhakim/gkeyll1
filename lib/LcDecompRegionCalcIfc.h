/**
 * @file	LcDecompRegionCalcIfc.h
 *
 * @brief	Base class for decomposition algorithms.
 */

#ifndef LC_DECOMP_REGION_CALC_IFC_H
#define LC_DECOMP_REGION_CALC_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegion.h>

namespace Lucee
{
/**
 * This class defines interface for different decomposition
 * algorithms. The actual algorithms are implemented by derived
 * classes and take the global region and the number of sub-regions as
 * input and compute the needed decomposition.
 */
  template <unsigned NDIM>
  class DecompRegionCalcIfc
  {
    public:
/**
 * Calculate decomposition adding subregions into decompRgn object.
 *
 * @param nrgns Number of regions to create.
 * @param decompRgn On output this contains decomposition.
 */
      void calcDecomp(unsigned nrgns, Lucee::DecompRegion<NDIM>& decompRgn);

    protected:

/**
 * Method used to compute actual decomposition. Derived classes should
 * provide this method and use `addRegion` method to add computed
 * regions.
 *
 * @param nrgns Number of regions to create.
 * @param globalRgn Global region to decompose.
 */
      virtual void decompose(unsigned nrgs,
        const Lucee::Region<NDIM, int>& globalRgn) = 0;

/**
 * Append sub-region to decomposition. This method should be called by
 * derived classes to add sub-regions.
 *
 * @param subRgn Region to add.
 */
      void addRegion(const Lucee::Region<NDIM, int>& subRgn);

    private:
/** Pointer to decomp region box: this is so addRegion can access it */
      Lucee::DecompRegion<NDIM> *decompRgnPtr;
  };
}

#endif // LC_DECOMP_REGION_CALC_IFC_H

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
#include <LcBasicObj.h>
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
  class DecompRegionCalcIfc : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new decomposition object */
      DecompRegionCalcIfc();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Set specified direction as periodic.
 *
 * @param dir Direction to set as periodic.
 * @param isp Flag to indicate if direction is periodic.
 */
      void setPeriodicDir(unsigned dir, bool isp = true);

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
/** Periodic directions */
      bool isPeriodic[NDIM];

/**
 * Account for periodic directions by duplicating boxes appropriately.
 */
      void handlePeriodicDirs();
  };
}

#endif // LC_DECOMP_REGION_CALC_IFC_H

/**
 * @file	LcRectYeeInterpolationUpdater.h
 *
 * @brief	Compute 2nd order central-differences on a rectangular grid.
 */

#ifndef LC_RECT_YEE_INTERPOLATION_UPDATER_H
#define LC_RECT_YEE_INTERPOLATION_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to compute second-order central difference of a field
 * defined on a rectangular grid.
 */
  template <unsigned NDIM>
  class RectYeeInterpolationUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new central-difference updater.
 */
      RectYeeInterpolationUpdater();
      
/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of updater.
 */
      Lucee::UpdaterStatus update(double t);

/**
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass. Inside the implementation of this method the
 * derived class must make a sequence of appendInpVarType() and
 * appendOutVarType() calls to declare the input/output data structure
 * types.
 */
      void declareTypes();

    private:
/**
 * Compute central difference of field.
 *
 * @param inFld Field to compute CD of.
 * @param cdFld Difference output.
 */
      void computeCentralDifference1D(const Lucee::Field<NDIM, double>& inFld, 
        Lucee::Field<NDIM, double>& cdFld);

/**
 * Compute central difference of field.
 *
 * @param inFld Field to compute CD of.
 * @param cdFld Difference output.
 */
      void computeCentralDifference2D(const Lucee::Field<NDIM, double>& inFld, 
        Lucee::Field<NDIM, double>& cdFld);

/**
 * Compute central difference of field.
 *
 * @param inFld Field to compute CD of.
 * @param cdFld Difference output.
 */
      void computeCentralDifference3D(const Lucee::Field<NDIM, double>& inFld, 
        Lucee::Field<NDIM, double>& cdFld);

      /** to interpolate to centre or edge */
      bool fwd; 
      /** Whether a Dual Yee cell or regular Yee cell is being used */
      bool dual; 
      /** Extra cells to update outside of interior */
      int ghostUpdates[2];
  };
}

#endif // LC_RECT_SECOND_ORDER_CENTRAL_DIFF_UPDATER_H

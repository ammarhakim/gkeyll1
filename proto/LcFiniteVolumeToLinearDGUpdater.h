/**
 * @file	LcFiniteVolumeToLinearDGUpdater.h
 *
 * @brief	Updater to take FV cell-center values and put them on a DG field
 */

#ifndef LC_FINITE_VOLUME_TO_LINEAR_DG_UPDATER_H
#define LC_FINITE_VOLUME_TO_LINEAR_DG_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDynVector.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to copy FV cell-center values and put them on a linear DG
 * field.
 */
  template <unsigned NDIM>
  class FiniteVolumeToLinearDGUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      FiniteVolumeToLinearDGUpdater();

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
/** Component to copy */
      unsigned component;
/** Flag to indicate if end node are to be extrapolated */
      bool extrapolateNodes;
  };
}

#endif // LC_FINITE_VOLUME_TO_LINEAR_DG_UPDATER_H

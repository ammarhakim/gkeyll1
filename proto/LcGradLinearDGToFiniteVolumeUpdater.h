/**
 * @file	LcGradLinearDGToFiniteVolumeUpdater.h
 *
 * @brief	Updater to take FV cell-center values and put them on a DG field
 */

#ifndef LC_GRAD_LINEAR_DG_TO_FINITE_VOLUME_UPDATER_H
#define LC_GRAD_LINEAR_DG_TO_FINITE_VOLUME_UPDATER_H

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
  class GradLinearDGToFiniteVolumeUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      GradLinearDGToFiniteVolumeUpdater();

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
  };
}

#endif // LC_GRAD_LINEAR_DG_TO_FINITE_VOLUME_UPDATER_H

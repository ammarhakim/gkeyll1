/**
 * @file	LcThreeWaveInteractSrcUpdater.h
 *
 * @brief	Three wave interaction source updater.
 */

#ifndef LC_THREE_WAVE_INTERACT_SRC_UPDATER_H
#define LC_THREE_WAVE_INTERACT_SRC_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Update source terms needed in three-wave interaction problem. This
 * system is give by the coupled equations
 *
 * a_t + a_z = -b f
 * b_t - b_z = a f*
 * f_t = -a b*
 *
 * Only source terms are updated in this updater.
 */
  class ThreeWaveInteractSrcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

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
  };
}

#endif // LC_THREE_WAVE_INTERACT_SRC_UPDATER_H

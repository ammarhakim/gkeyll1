/**
 * @file	LcTwoFluidMomentumRelaxSrcUpdater.h
 *
 * @brief	Updater to apply momentum relaxation from inter-species collisions
 */

#ifndef LC_TWO_FLUID_MOMENTUM_RELAX_SRC_UPDATER_H
#define LC_TWO_FLUID_MOMENTUM_RELAX_SRC_UPDATER_H

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
 * Update fluid momentum from relaxation terms due to inter-species
 * collisions.
 */
  template <unsigned NDIM>
  class TwoFluidMomentumRelaxSrcUpdater : public Lucee::UpdaterIfc
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
/** Electron collision frequency */
      double elcNu;
  };
}

#endif // LC_TWO_FLUID_MOMENTUM_RELAX_SRC_UPDATER_H

/**
 * @file	LcFiveMomentNumDensityRelax.h
 *
 * @brief	Updater to apply momentum relaxation from inter-species collisions
 */

#ifndef LC_FIVE_MOMENT_NUM_DENSITY_RELAX_H
#define LC_FIVE_MOMENT_NUM_DENSITY_RELAX_H

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
 * Relax number density to enforce quasineutrality in two-fluid
 * equations.
 */
  template <unsigned NDIM>
  class FiveMomentNumDensityRelax : public Lucee::UpdaterIfc
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
/** Electron mass */
      double elcMass;
/** Ion mass */
      double ionMass;
/** Gas gamma */
      double gasGamma;
  };
}

#endif // LC_FIVE_MOMENT_NUM_DENSITY_RELAX_H

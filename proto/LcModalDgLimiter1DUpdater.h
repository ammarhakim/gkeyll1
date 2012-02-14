/**
 * @file	LcModalDgLimiter1DUpdater.h
 *
 * @brief	Updater to apply limiters to DG solution.
 */

#ifndef LC_MODAL_DG_LIMITER_1D_UPDATER_H
#define LC_MODAL_DG_LIMITER_1D_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to apply limiter to solution computed using modal DG
 * scheme.
 */
  class ModalDgLimiter1DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG limiter 1D */
      ModalDgLimiter1DUpdater();

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
/** Equation to solve */
      Lucee::HyperEquation *equation;
/** Number of equations */
      unsigned meqn;
/** Number of basis functions to use */
      double numBasis;
/** Factor for slope correction */
      double Mfact;

/**
 * Modified min-mod function of three variables. Usually minmod
 * function does care about the order of the parameters. In this
 * function, however, it does: the first parameter is the slope
 * returned if it lies in some interval.
 *
 * @param a First parameter.
 * @param b Second parameter.
 * @param c Third parameter.
 * @param dx Cell spacing
 */
      double modifiedMinMod(double a, double b, double c, double dx) const;
  };
}

#endif // LC_MODAL_DG_1D_UPDATER_H

/**
 * @file	LcMusclHancock1DUpdater.h
 *
 * @brief	Solver for 1D Euler equations using MUSCl-Hancock scheme.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_MUSCL_HANCOCK_1D_UPDATER_H
#define LC_MUSCL_HANCOCK_1D_UPDATER_H

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
 * Updater to solve 1D Euler equations using MUSCL-Hancock scheme.
 */
  class MusclHancock1DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new combiner.
 */
      MusclHancock1DUpdater();

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
/** Limiter to use */
      unsigned limiter;
/** Gas adibatic constant */
      double gas_gamma;
/** CFL number */
      double cfl;
/** Field to store slopes */
      Lucee::Field<1, double> slopes;
/** Field to store predicted variables */
      Lucee::Field<1, double> predict;

/**
 * Averaging function. This returns a "limited" average that may
 * prevent non-physical oscillations depending on the specified
 * limiter to use.
 *
 * @param a First value in average.
 * @param b Second value in average.
 * @return average, possibly limited.
 */
      double limaverage(double a, double b);

/**
 * Compute primitive variables from conserved variables.
 *
 * @param cv Conserved variables.
 * @param pv (out) Primitive variables.
 */
      void calcPrimVars(const double *cv, double *pv);
      
  };
}

#endif // LC_MUSCL_HANCOCK_1D_UPDATER_H

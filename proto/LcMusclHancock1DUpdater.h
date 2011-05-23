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

/**
 * Averaging function. This returns a "limited" average that may
 * prevent non-physical oscillations depending on the specified
 * limiter to use.
 */
      double limave(double a, double b);
  };
}

#endif // LC_MUSCL_HANCOCK_1D_UPDATER_H

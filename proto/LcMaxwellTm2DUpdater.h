/**
 * @file	LcMaxwellTm2DUpdater.h
 *
 * @brief	Solver for transverse-magnetic Maxwell equations in 2D.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_MAXWELL_TM_2D_UPDATER_H
#define LC_MAXWELL_TM_2D_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to solve transverse-magnetic Maxwell equations in 2D. This
 * is a test updater and is not meant for production simulations.
 */
  class MaxwellTm2DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new combiner.
 */
      MaxwellTm2DUpdater();

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
/** Speed of light */
      double c0;
/** Correction speed factor */
      double gamma;
  };
}

#endif // LC_MAXWELL_TM_2D_UPDATER_H

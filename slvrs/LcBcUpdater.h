/**
 * @file	LcBcUpdater.h
 *
 * @brief	Apply boundary conditions.
 */

#ifndef LC_BC_UPDATER_H
#define LC_BC_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to apply boundary conditions on structured grids. This
 * updater takes a list of BCs to apply and applies them to the
 * supplied vectors.
 */
  template <unsigned NDIM>
  class BcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new curl updater.
 */
      BcUpdater();

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
/** Direction to apply boundary condtion */
      unsigned dir;
/** Edge to apply boundary condition */
      unsigned edge;
  };
}

#endif // LC_BC_UPDATER_H

/**
 * @file	LcGridOdePointIntegratorUpdater.h
 *
 * @brief	Updater to integrate ODEs on a grid.
 */

#ifndef LC_GRID_ODE_POINT_UPDATER_H
#define LC_GRID_ODE_POINT_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridOdePointIntegrator.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to integrate ODEs on a grid.
 */
  template <unsigned NDIM>
  class GridOdePointIntegratorUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new integrator.
 */
      GridOdePointIntegratorUpdater();

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
/** Ode integrator to use */
      Lucee::GridOdePointIntegrator<NDIM> *integrator;
  };
}

#endif // LC_GRID_ODE_POINT_UPDATER_H

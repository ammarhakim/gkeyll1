/**
 * @file	LcStairSteppedBcUpdater.h
 *
 * @brief	Base class for boundary conditions for stair-stepped boundaries.
 */

#ifndef LC_STAIR_STEPPED_BC_UPDATER_H
#define LC_STAIR_STEPPED_BC_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Base class for boundary conditions for stair-stepped
 * boundaries.
 */
  template <unsigned NDIM>
  class StairSteppedBcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      StairSteppedBcUpdater();

/**
 * Delete updater.
 */
      virtual ~StairSteppedBcUpdater();

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
/** Direction to apply boundary conditon */
      unsigned dir;
/** Boundary conditions to apply */
      std::vector<Lucee::BoundaryCondition*> bcList;
/** Pointer to in/out field */
      Lucee::Field<NDIM, double> *inOut;
/** Field to store information about boundary */
      Lucee::Field<NDIM, double> *ssBnd;
  };
}

#endif // LC_STAIR_STEPPED_BC_UPDATER_H

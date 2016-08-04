/**
 * @file	LcMultplyFieldsUpdater.h
 *
 * @brief	Updater to multiply an arbitrary number of fields to each other
 */

#ifndef LC_MULTIPLY_FIELDS_UPDATER_H
#define LC_MULTIPLY_FIELDS_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDynVector.h>
#include <LcField.h>
#include <LcUpdaterIfc.h>
#include <LcNodalFiniteElementIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to integrate nodal DG/CG field over complete domain.
 */
  template <unsigned NDIM>
  class MultiplyFieldsUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      MultiplyFieldsUpdater();

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

#endif // LC_MULTIPLY_FIELDS_UPDATER_H

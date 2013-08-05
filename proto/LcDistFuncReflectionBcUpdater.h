/**
 * @file	LcDistFuncReflectionBcUpdater.h
 *
 * @brief	Applies particle refection BCs to distribution function
 */

#ifndef LC_DIST_FUNC_REFLECTION_BS_UPDATER_H
#define LC_DIST_FUNC_REFLECTION_BS_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class DistFuncReflectionBcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      DistFuncReflectionBcUpdater();

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
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis;
/** Edges to apply boundary condition */
      unsigned edge;
  };
}

#endif // LC_DIST_FUNC_REFLECTION_BS_UPDATER_H

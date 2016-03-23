/**
 * @file	LcSOLSetPotentialAtBoundary.h
 *
 * @brief	Sets Dirchlet boundary conditions on potential for 5D SOL simulations
 */

#ifndef LC_SOL_SET_POTENTIAL_AT_BOUNDARY_H
#define LC_SOL_SET_POTENTIAL_AT_BOUNDARY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcDynVector.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater
 */
  class SOLSetPotentialAtBoundary : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new solver */
      SOLSetPotentialAtBoundary();
/** Destructor */
      virtual ~SOLSetPotentialAtBoundary();

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
/** Pointer to 3D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
  };
}

#endif // LC_SOL_SET_POTENTIAL_AT_BOUNDARY_H

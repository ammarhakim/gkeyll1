/**
 * @file	LcDistFuncReflectionBcUpdater.h
 *
 * @brief	Applies particle refection BCs to (electron) distribution function
 */

#ifndef LC_REFLECTION_BOUNDARY_CONDITION_H
#define LC_REFLECTION_BOUNDARY_CONDITION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class ReflectionBoundaryCondition : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      ReflectionBoundaryCondition();

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
      Lucee::NodalFiniteElementIfc<3> *nodalBasis;
/** Flag to indicate if BCs should be applied to left edge */
      bool applyLeftEdge;
/** Flag to indicate if BCs should be applied to right edge */
      bool applyRightEdge;
/** Mapping for 180 degree rotations */
      std::vector<unsigned> rotMap;
  };
}

#endif // LC_REFLECTION_BOUNDARY_CONDITION_H

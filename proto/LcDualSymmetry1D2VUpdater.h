/**
 * @file	LcDualSymmetry1D2VUpdater.h
 *
 * @brief	Simple updater that takes in a function W(y,py,px) and returns
 * and updater with the result W' = 0.5*[W(y,py,px) + W(y,-py,px)]
 */

#ifndef LC_DUAL_SYMMETRY_1D_2V_UPDATER_H
#define LC_DUAL_SYMMETRY_1D_2V_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcStructGridField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to evaluate the diffusion term in the Lenard-Bernstein
 * collision operator.
 */
  class DualSymmetry1D2VUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      DualSymmetry1D2VUpdater();

/**
 * Destroy updater.
 */
      ~DualSymmetry1D2VUpdater();

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
/** Mapping for 180 degree rotations */
      std::vector<unsigned> rotMap;
  };
}

#endif // LC_DUAL_SYMMETRY_1D_2V_UPDATER_H

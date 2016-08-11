/**
 * @file	LcSimpleSmoothToC0Updater.h
 *
 * @brief	Updater to smooth a field by averaging nodal values on shared edges
 */

#ifndef LC_SIMPLE_SMOOTH_TO_C0_UPDATER_H
#define LC_SIMPLE_SMOOTH_TO_C0_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
  template <unsigned NDIM>
  class SimpleSmoothToC0Updater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG updater */
      SimpleSmoothToC0Updater();

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
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;
/**
 * Structure to store node numbers on edges.
 */
      struct EdgeNodeNums
      {
/** Node numbers */
          std::vector<int> nums;
      };

/** Vector to store lower node numbers */
      EdgeNodeNums lowerNodeNums[NDIM];
/** Vector to store upper node numbers */
      EdgeNodeNums upperNodeNums[NDIM];
/** Temporary flag to keep track of what polynomial order of element we are smoothing */
      int polyOrder;
  };
}

#endif // LC_SIMPLE_SMOOTH_TO_C0_UPDATER_H

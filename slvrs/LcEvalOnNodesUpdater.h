/**
 * @file	LcEvalOnNodesUpdater.h
 *
 * @brief	Evaluate function on nodal field.
 */

#ifndef LC_EVAL_ON_NODES_UPDATER_H
#define LC_EVAL_ON_NODES_UPDATER_H

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
 * Updater to set a nodal FE field from supplied Lua function
 */
  template <unsigned NDIM>
  class EvalOnNodesUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      EvalOnNodesUpdater();

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
/** Reference to function to project */
      int fnRef;
/** Flag to indicate if nodes are shared */
      bool sharedNodes;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;

/**
 * Evaluate function at specified location and fill output array with
 * result.
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param nc Matrix with nodal cooridates.
 * @param nn Number number
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        const Lucee::Matrix<double> nc, double nn, std::vector<double>& res);
  };
}

#endif // LC_EVAL_ON_NODES_UPDATER_H

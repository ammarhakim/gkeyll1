/**
 * @file	LcConstructLinearOperatorMatrix.h
 *
 * @brief	Given df/dt and an index, puts results into a column of a sparse matrix.
 * This matrix, when complete, is the linear matrix L such that df/dt = L*f.
 */

#ifndef LC_CONSTRUCT_LINEAR_OPERATOR_MATRIX_H
#define LC_CONSTRUCT_LINEAR_OPERATOR_MATRIX_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

// eigen includes
#include <Eigen/SparseCore>

namespace Lucee
{
/**
 * Updater to set a nodal FE field from supplied Lua function
 */
  template <unsigned NDIM>
  class ConstructLinearOperatorMatrix : public Lucee::UpdaterIfc
  {
// Number of components for coordinate arrays etc.
      static const unsigned NC = Lucee::CDIM<NDIM>::N;

    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      ConstructLinearOperatorMatrix();

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
/** Reference to function indicating write location */
      int fnRef;
/** Stores total number of degrees of freedom to store per column */
      int totalNodes;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;
/** Stores list of entries to insert into sparse matrix */
      std::vector<Eigen::Triplet<double> > tripletList;
/**
 * Evaluate function at specified location and fill output array with
 * result.
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        std::vector<double>& res);
  };
}

#endif // LC_CONSTRUCT_LINEAR_OPERATOR_MATRIX_H

/**
 * @file	LcNodalPositiveFilterUpdater.h
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 */

#ifndef LC_NODAL_POSITIVE_FILTER_UPDATER_H
#define LC_NODAL_POSITIVE_FILTER_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to apply positivity filer to DG solution.
 */
  template <unsigned NDIM>
  class NodalPositiveFilterUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      NodalPositiveFilterUpdater();

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

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      //static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    private:
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;
/** Number of local nodes */
      unsigned nlocal;
/** Number of equations */
      unsigned meqn;
/** Equation to solve */
      Lucee::HyperEquation *equation;
/** Weights for quadrature */
      std::vector<double> weights;

/** Enum for operation type */
      enum Operation { OP_FLATTEN, OP_FILTER };
/** Operation to perform */
      Operation opType;

/**
 * Compute cell averages from nodal values.
 *
 * @param [in] qIn Input nodal values.
 * @param [out] qAvg On output, average values.
 */
      void calcAverage(const Lucee::ConstFieldPtr<double>& qIn, std::vector<double>& qAvg);

/**
 * Reset nodal values to specified values.
 *
 * @param q On output, all nodal values in this will be replace.
 * @param qAvg Input values to reset to.
 */
      void resetAllNodes(Lucee::FieldPtr<double>& q, const std::vector<double>& qAvg);
  };
}

#endif // LC_NODAL_POSITIVE_FILTER_UPDATER_H

/**
 * @file	LcSOLLenardBernsteinScaleCell5DUpdater.h
 *
 * @brief	Accumulates correct amount of an input distribution function (assumed to be a diffusion term)
 * to an existing distribution so that the desired amount of energy in each cell is achieved
 */

#ifndef LC_SOL_LENARD_BERNSTEIN_SCALE_CELL_5D_UPDATER_H
#define LC_SOL_LENARD_BERNSTEIN_SCALE_CELL_5D_UPDATER_H

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

// std includes
#include <vector>

// eigen includes
#include <Eigen/Core>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class SOLLenardBernsteinScaleCell5DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOLLenardBernsteinScaleCell5DUpdater();

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
/** Pointer to phase space basis functions to use */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to configuration space basis functions */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** CFL number */
      double cfl;
/** Reference to function specifying alpha */
      int fnRef;
/** mass of species */
      double speciesMass;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
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

#endif // LC_SOL_LENARD_BERNSTEIN_SCALE_CELL_5D_UPDATER_H

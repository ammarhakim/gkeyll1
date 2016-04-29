/**
 * @file	LcSOLTotalEnergyCalc.h
 *
 * @brief	Diagnostic to calculate total energy in 5D SOL simulation
 */

#ifndef LC_SOL_TOTAL_ENERGY_CALC_H
#define LC_SOL_TOTAL_ENERGY_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcDynVector.h>

// eigen includes
#include <Eigen/Core>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class SOLTotalEnergyCalc : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOLTotalEnergyCalc();

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
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/**
 * Volume interpolation matrix in 5d
 */
      Eigen::MatrixXd volQuad5d;
/**
 * Quadrature weights for points on the upper z surface
 */
      std::vector<double> volWeights5d;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_SOL_TOTAL_ENERGY_CALC_H

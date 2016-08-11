/**
 * @file	LcSOLUpperXPotentialBcUpdater.h
 *
 * @brief	
 */

#ifndef LC_SOL_UPPER_X_POTENTIAL_BC_UPDATER_H
#define LC_SOL_UPPER_X_POTENTIAL_BC_UPDATER_H

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
 * Updater to compute zonal average of phi(x,y,z)
 */
  class SOLUpperXPotentialBcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new solver */
      SOLUpperXPotentialBcUpdater();
/** Destructor */
      virtual ~SOLUpperXPotentialBcUpdater();

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
/** Pointer to 2D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** When multiplied with vector of solution in a cell, gives y-integrated value on upper x surface */
      Eigen::VectorXd integrationMatrixUpper;
/** When multiplied with vector of solution in a cell, gives y-integrated value on lower x surface */
      Eigen::VectorXd integrationMatrixLower;
/** Keeps track of nodes to write data to on upper surface in x */
      std::vector<int> upperEdgeNodeNums;
/** Keeps track of nodes to write data to on lower surface in x */
      std::vector<int> lowerEdgeNodeNums;
/** Flag indicating whether or not to apply the same bc to lower x edge */
      bool applyLowerEdge;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_SOL_UPPER_X_POTENTIAL_BC_UPDATER_H

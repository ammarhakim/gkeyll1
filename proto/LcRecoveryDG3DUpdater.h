/**
 * @file	LcRecoveryDG3DUpdater.h
 *
 * @brief	Updater to perform recovery DG calculation for 3d problems.
 * Currently supports calculation of second and fourth derivatives of polyOrder = 1 input field
 */

#ifndef LC_RECOVERY_DG_3D_UPDATER_H
#define LC_RECOVERY_DG_3D_UPDATER_H

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

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to evaluate the diffusion term in the Lenard-Bernstein
 * collision operator.
 */
  class RecoveryDG3DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      RecoveryDG3DUpdater();

/**
 * Destroy updater.
 */
      ~RecoveryDG3DUpdater();

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
/** CFL number */
      double cfl;
/** Flag to indicate if to compute only increments */
      bool onlyIncrement;
/** if false, we will ignore cfl checks */
      bool checkTimeStepSize;
/** Basis function order */
      int polyOrder;
/** Order of derivative to calculate */
      int derivOrder;
/** Factor in front of derivative */
      double alpha;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis;
/** Eigen matrices for recovery calculation */
      Eigen::MatrixXd lowerMat;
      Eigen::MatrixXd upperMat;
      Eigen::MatrixXd selfCenter;
      Eigen::MatrixXd lowerCenter;
      Eigen::MatrixXd upperCenter;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_RECOVERY_DG_3D_UPDATER_H

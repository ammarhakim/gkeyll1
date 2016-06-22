/**
 * @file	LcZonalVelocity1DUpdater.h
 *
 * @brief	
 */

#ifndef LC_ZONAL_VELOCITY_1D_UPDATER_H
#define LC_ZONAL_VELOCITY_1D_UPDATER_H

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
  class ZonalVelocity1DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      ZonalVelocity1DUpdater();

/**
 * Destroy updater.
 */
      ~ZonalVelocity1DUpdater();

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
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis;
/** Grad stiffness matrix */
      Eigen::MatrixXd gradStiffnessMatrix;
      Eigen::MatrixXd lowerSurfInterpMatrix;
      Eigen::MatrixXd upperSurfInterpMatrix;
      Eigen::MatrixXd massMatrixInv;
      std::vector<double> lowerSurfWeights;
      std::vector<double> upperSurfWeights;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_RECOVERY_DG_1D_UPDATER_H

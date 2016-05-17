/**
 * @file	LcLenardBernsteinDiff5DUpdater.h
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator for 3D2V problems.
 */

#ifndef LC_LENARD_BERNSTEIN_DIFF_5D_UPDATER_H
#define LC_LENARD_BERNSTEIN_DIFF_5D_UPDATER_H

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
  class LenardBernsteinDiff5DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      LenardBernsteinDiff5DUpdater();

/**
 * Destroy updater.
 */
      ~LenardBernsteinDiff5DUpdater();

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
/** Diffusion coefficient */
      double alpha;
/** CFL number */
      double cfl;
/** Flag to indicate if to compute only increments */
      bool onlyIncrement;
/** Reference to function specifying alpha */
      int fnRef;
/** species mass for vPerp to mu transformation */
      double speciesMass;
/** Pointer to 5d basis functions */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to 3d basis functions */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** List of matrices on each lower face */
      std::vector<Eigen::MatrixXd > lowerMat;
/** List of matrices on each upper face */
      std::vector<Eigen::MatrixXd > upperMat;
/** Precomputed matrices accounting for how mu is a radial coordinate */
      Eigen::MatrixXd iMatDiffusionTimesMu;
      Eigen::MatrixXd lowerMatDiffusionTimesMu;
      Eigen::MatrixXd upperMatDiffusionTimesMu;
/** Vectors to store eigen matrices */
      std::vector<Eigen::MatrixXd> upperCenter;
      std::vector<Eigen::MatrixXd> selfCenter;
      std::vector<Eigen::MatrixXd> lowerCenter;
/** Additional matrices needed for diffusion in mu coordinate */
      Eigen::MatrixXd upperCenterTimesMu;
      Eigen::MatrixXd selfCenterTimesMu;
      Eigen::MatrixXd lowerCenterTimesMu;
/** Weights for 3D volume gaussian quadrature points */
      std::vector<double> gaussVolWeights;
/** 
 * 3D interpolation matrix for bringing quantities from nodal locations to volume
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpVolMatrix;
/**
 * When multiplied by the solution in a 3d cell, returns the integrated cell quantity
 * To get average, one should divide by the volume of the cell
 */
      Eigen::VectorXd mom0Vector;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::Matrix<double, 32, 32>& destinationMatrix);
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

#endif // LC_LENARD_BERNSTEIN_DIFF_5D_UPDATER_H

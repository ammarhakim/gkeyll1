/**
 * @file	LcSOLElectronDensityInitialization.h
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 */

#ifndef LC_SOL_ELECTRON_DENSITY_INITIALIZATION_H
#define LC_SOL_ELECTRON_DENSITY_INITIALIZATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve hyperbolic equations using a nodal discontinous
 * Galerkin scheme.
 */
  class SOLElectronDensityInitialization : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      SOLElectronDensityInitialization();

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
      Lucee::NodalFiniteElementIfc<1> *nodalBasis;
/** Value of k_perp0*rho_s */
      double kPerpTimesRho;
/** Weights for gaussian quadrature points */
      std::vector<double> gaussWeights;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpMatrix;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);

      // Testing GSL functions
      struct fParams {
        double averageIonDensity;
        double ionDensityAtPoint;
      };
      double computeF(double x, void *params);
      double computeFPrime(double x, void *params);
      void computeFAndFPrime(double x, void *params, double *y, double *dy);
      
      static double computeF_Wrapper(double x, void *params);
      static double computeFPrime_Wrapper(double x, void *params);
      static void computeFAndFPrime_Wrapper(double x, void *params, double *y, double *dy);
  };
}

#endif // LC_SOL_ELECTRON_DENSITY_INITIALIZATION_H

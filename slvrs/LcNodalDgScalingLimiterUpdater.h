/**
 * @file	LcNodalDgScalingLimiterUpdater.h
 *
 * @brief	Updater to enforce positivity preservation of the distribution function in Vlasov simulations
 */

#ifndef LC_VLASOV_POSITIVITY_UPDATER_H
#define LC_VLASOV_POSITIVITY_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcCDIM.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to compute moments of distribution function f(x,v)
 */
  template <unsigned CDIM, unsigned VDIM>
  class NodalDgScalingLimiterUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned PNC = Lucee::CDIM<CDIM+VDIM>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned CNC = Lucee::CDIM<CDIM>::N;

/** Create new modal DG solver in CDIM dimensions */
      NodalDgScalingLimiterUpdater();

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
/** Pointer to nodal phase space (CDIM+VDIM) basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM+VDIM> *phaseBasis;
/** Pointer to nodal configuration space (CDIM) basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;
/** Vector of configuration space quadrature weights to compute average */
      Eigen::VectorXd volWeights;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_VLASOV_POSITIVITY_UPDATER_H

/**
 * @file	LcModalDg1DDiffusionUpdater.h
 *
 * @brief	Updater to solver 1D hyper diffusion equations using modal DG scheme
 */

#ifndef LC_MODAL_DG_1D_HYPER_DIFFUSION_UPDATER_H
#define LC_MODAL_DG_1D_HYPER_DIFFUSION_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

// etc includes
#include <quadrule.hpp>

namespace Lucee
{
/**
 * Updater to solve 1D hyperbolic equations using modal DG
 * scheme. This updater computes a first-order Euler update
 *
 * qNew = q + dt * L(q)
 *
 * where dt is the specified time-step and L(q) the RHS of the
 * semi-discrete equation. Using the output of this updater any RK
 * scheme can be easily implemented in Lua itself.
 */
  class ModalDg1DHyperDiffusionUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      ModalDg1DHyperDiffusionUpdater();

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
/** CFL number to use */
      double cfl;
/** Maximum CFL number */
      double cflm;
/** Number of basis functions to use */
      unsigned numBasis;
/** Diffusion coefficient */
      double diffCoef;
/** Normalization values */
      std::vector<double> normCoeff;
/** Inverse of recovery matrix */
      Eigen::MatrixXd recoveryMatrixInv;
/** Precomputed volume integrals */
      Eigen::MatrixXd recoveryVolume;
/** 2nd Derivative of Legendre Polynomial evaluated at edges -1 and +1
  * Col 0 = -1 edge, Col 1 = +1 edge (then normalied to physical space)
  */
      Eigen::MatrixXd secondDerivEdgeEvals;
/** 3rd Derivative of Legendre Polynomial evaluated at edges -1 and +1
  * Col 0 = -1 edge, Col 1 = +1 edge (then normalied to physical space)
  */
      Eigen::MatrixXd thirdDerivEdgeEvals;
/** Flag to indicate if only increments should be computed */
      bool onlyIncrement;
/**
 * Evaluate the second derivative of a legendre polynomial (order 'order') at
 * the location x
 */
      double legendrePoly2ndDeriv(int order, double x);
/**
 * Returns the coefficients of "1" and "x" terms of recovery polynomial.
 * col 1 = f0, col 2 = f1, each row is f's for one eqn
 */
      void computeRecoveryPolynomial(const Lucee::ConstFieldPtr<double>& qLeft,
        Lucee::ConstFieldPtr<double>& qRight, std::vector<double>& boundaryTerms);
  };
}

#endif // LC_MODAL_DG_1D_HYPER_DIFFUSION_UPDATER_H

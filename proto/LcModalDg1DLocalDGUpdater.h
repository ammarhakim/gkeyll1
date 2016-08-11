/**
 * @file	LcModalDg1DLocalDGUpdater.h
 *
 * @brief	Updater to solver 1D diffusion equations using the local DG scheme
 */

#ifndef LC_MODAL_DG_1D_LOCAL_DG_UPDATER_H
#define LC_MODAL_DG_1D_LOCAL_DG_UPDATER_H

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
  class ModalDg1DLocalDGUpdater : public Lucee::UpdaterIfc
  {
/** Type of flux pairs to use */
    enum FluxPair { RIGHT, LEFT, AVERAGE };
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      ModalDg1DLocalDGUpdater();

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
/** Sqrt of diffusion coefficient */
      double diffCoefSqrt;
/** Normalization values */
      std::vector<double> normCoeff;
/** Precomputed volume integrals */
      Eigen::MatrixXd volumeMatrix;
/** Flag to indicate if only increments should be computed */
      bool onlyIncrement;
/** Type of flux pair to use */
      FluxPair fluxPair;
  };
}

#endif // LC_MODAL_DG_1D_LOCAL_DG_UPDATER_H

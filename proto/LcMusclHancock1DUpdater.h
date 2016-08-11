/**
 * @file	LcMusclHancock1DUpdater.h
 *
 * @brief	Solver for 1D Euler equations using MUSCl-Hancock scheme.
 */

#ifndef LC_MUSCL_HANCOCK_1D_UPDATER_H
#define LC_MUSCL_HANCOCK_1D_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to solve 1D Euler equations using MUSCL-Hancock scheme.
 */
  class MusclHancock1DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new combiner.
 */
      MusclHancock1DUpdater();

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
/** Limiter to use */
      unsigned limiter;
/** Gas adibatic constant */
      double gas_gamma;
/** CFL number */
      double cfl;
/** Factor for use in EPSILON limiter */
      double epsFac;
/** Field to store slopes */
      Lucee::Field<1, double> slopes;
/** Field to store predicted variables */
      Lucee::Field<1, double> predict;
/** Field to store primitive variables */
      Lucee::Field<1, double> prim;

/**
 * Averaging function. This returns a "limited" average that may
 * prevent non-physical oscillations depending on the specified
 * limiter to use.
 *
 * @param a First value in average.
 * @param b Second value in average.
 * @return average, possibly limited.
 */
      double limaverage(double a, double b);

/**
 * Compute primitive variables from conserved variables.
 *
 * @param cv Conserved variables.
 * @param pv (out) Primitive variables.
 */
      void calcPrimVars(const Lucee::Field<1, double>& cv, Lucee::Field<1, double> &pv);

/**
 * Calculate numerical fluxes given primitive variables on left and
 * right of edge.
 *
 * @param pvl Primitive variables on left of edge.
 * @param pvr Primitive variables on right of edge.
 * @param nf (out) Numerical flux.
 */
      void calcNumericalFlux(const std::vector<double>& pvl, const std::vector<double> &pvr,
        std::vector<double>& nf);

/**
 * Compute fluxes given primitive variables.
 *
 * @param pv Primitive variables.
 * @param flux (out) Flux.
 */
      void calcFlux(const std::vector<double>& pv, std::vector<double>& flux);

/**
 * Compute conserved variables given primitive variables.
 *
 * @param pv Primitive variables.
 * @param cv (out) Conserved variables.
 */
      void calcConsVars(const std::vector<double>& pv, std::vector<double>& cv);

/**
 * Project vector on left eigenvectors of flux Jacobian.
 *
 * @param pv Primitive variables with which to compute flux Jacobian.
 * @param vec Vector to project.
 * @param coeff (out) Vector of coefficients.
 */
      void projectOnLeftEigenvectors(const double *pv, const double *vec, double *coeff);

/**
 * Reconstruct vector by weighted sum of right eigenvectors.
 *
 * @param pv Primitive variables with which to compute flux Jacobian.
 * @param coeff Vector of coefficients.
 * @param vec (out) Output vector.
 */
      void reconWithRightEigenvectors(const double *pv, const double *coeff, double *vec);
  };
}

#endif // LC_MUSCL_HANCOCK_1D_UPDATER_H

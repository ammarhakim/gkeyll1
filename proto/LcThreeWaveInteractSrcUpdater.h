/**
 * @file	LcThreeWaveInteractSrcUpdater.h
 *
 * @brief	Three wave interaction source updater.
 */

#ifndef LC_THREE_WAVE_INTERACT_SRC_UPDATER_H
#define LC_THREE_WAVE_INTERACT_SRC_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>
#include <complex>

namespace Lucee
{
/** Integrator type */
  enum {TWI_RK4, TWI_IMPLICIT};

/**
 * Update source terms needed in three-wave interaction problem. This
 * system is give by the coupled equations
 *
 * a_t + a_z = -b f
 * b_t - b_z = a f*
 * f_t = -a b*
 *
 * Only source terms are updated in this updater.
 */
  class ThreeWaveInteractSrcUpdater : public Lucee::UpdaterIfc
  {
/** State of three-wave RHS */
      typedef std::vector<std::complex<double> > twstate;

    public:
/** Class id: this is used by registration system */
      static const char *id;

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

/**
 * Functor for use in boost::odeint integrators.
 *
 * @param x Input state.
 * @param dxdt RHS of ODE
 * @param t Time at which RHS is computed.
 */
      void operator() (const twstate &x, twstate &dxdt, const double t);

    private:
/** Relative tolerance for ODE solver */
      double relTol;
/** List of constants multiplying quadratic */
      std::complex<double> c[3];
/** Integrator type */
      int intType;

/** 
 * Integrator using iterative implicit scheme.
 *
 * @param dt Time-step to advance ODE by.
 * @param inp Initial condition.
 * @param out Updated solution.
 * @return Number of iterations taken.
 */
      unsigned stepImplicit(double dt, const std::vector<std::complex<double> >& inp,
        std::vector<std::complex<double> >& out);

/** 
 * Integrator using Runge-Kutta 4th order scheme.
 *
 * @param dt Time-step to advance ODE by.
 * @param inp Initial condition.
 * @param out Updated solution.
 * @return Number of iterations (sub-steps) taken.
 */
      unsigned stepRK4(double dt, const std::vector<std::complex<double> >& inp,
        std::vector<std::complex<double> >& out);

/**
 * Compare two numbers to specified tolerance (set in relTol)
 *
 * @param a First number to compare
 * @param b Second number to compare
 * @return true if within tolerance, false otherwise
 */
      bool epsCmp(double a, double b) const;
  };
}

#endif // LC_THREE_WAVE_INTERACT_SRC_UPDATER_H

/**
 * @file	LcPeriodicCollisionlessHeatFluxUpdater.h
 *
 * @brief	Updater to solve Poisson equations with periodic BCs
 */

#ifndef LC_PERIODIC_COLLISIONLESS_HEAT_FLUX_UPDATER_H
#define LC_PERIODIC_COLLISIONLESS_HEAT_FULX_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUpdaterIfc.h>

// FFTW includes
# include <fftw3.h>
#ifdef HAVE_MPI
#include <fftw3-mpi.h>
#endif

// std includes
#include <complex>
#include <vector>

namespace Lucee
{
/**
 * Updater to solve Poisson equations in 2D with periodic BCs.
 */
  template <unsigned NDIM>
  class PeriodicCollisionlessHeatFluxUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      PeriodicCollisionlessHeatFluxUpdater();

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
/** Copy of source array for use in FFTW. Also used to store output FFTW */
      std::vector<std::complex<double> > src_in_out;
/** Data for FFTW of solution. Also used to store output inverse FFTW */
      std::vector<std::complex<double> > sol_in_out;
/** X-direction wave-numbers */
      std::vector<double> kx;
/** Y-direction wave-numbers */
      std::vector<double> ky;
/** array to store abs(k) */
      std::vector<double> kabs;
/** Pointers for use in FFTW for source */
      fftw_complex *f_src_in_out;
/** Pointers for use in FFTW for solution */
      fftw_complex *f_sol_in_out;
/** FFTW-Plan for forward transform */
      fftw_plan f_plan;
/** FFTW-Plan for inverse transform */
      fftw_plan b_plan;
  };
}

#endif // LC_PERIODIC_POISSON_2D_UPDATER_H

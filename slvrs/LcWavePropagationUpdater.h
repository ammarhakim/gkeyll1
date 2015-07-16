/**
 * @file	LcWavePropagationUpdater.h
 *
 * @brief	Wave propagation solver.
 */

#ifndef LC_WAVE_PROPAGATION_UPDATER_H
#define LC_WAVE_PROPAGATION_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to solve hyperbolic equations using the wave-propagation
 * scheme.
 */
  template <unsigned NDIM>
  class WavePropagationUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Delete allocated memory */
      virtual ~WavePropagationUpdater();

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
/** Directions to update */
      std::vector<unsigned> updateDims;
/** Limiter to use */
      unsigned limiter;
/** Equation to solve */
      Lucee::HyperEquation *equation;
/** CFL number to use */
      double cfl;
/** Maximum CFL number */
      double cflm;
/** Flag to indicate if a lower in boundary a direction has flux
 * boundaries specified */
      bool hasLowerFluxBc[NDIM];
/** Flag to indicate if a upper in boundary a direction has flux
 * boundaries specified */
      bool hasUpperFluxBc[NDIM];
/** Fields to store positive fluctuations */
      std::vector<Lucee::Field<1, double>* > apdq;
/** Fields to store negative fluctuations */
      std::vector<Lucee::Field<1, double>* > amdq;
/** Fields to store speeds */
      std::vector<Lucee::Field<1, double>* > speeds;
/** Fields to store waves */
      std::vector<Lucee::Field<1, double>* > waves;
/** Fields to store second order corrections */
      std::vector<Lucee::Field<1, double>* > fs;
/** Flag to indicate if we there is an embedded boundary */
      bool hasSsBnd;
/** In/out field */
      Lucee::Field<NDIM, double> *inOut;
/** Flux boundary field */
      Lucee::Field<NDIM, double> *fluxBc;
/**
 * Apply limiters to waves.
 *
 * @param dir Direction in which limiter is being applied.
 * @param idx Index of cell (this is only partially valid: idx[dir] is not properly set)
 * @param waves [in/out] On input, waves to be limited. On output, limited waves.
 * @param speeds Wave speeds.
 */
      void applyLimiters(unsigned dir, int idx[NDIM],
        Lucee::Field<1, double>& waves, const Lucee::Field<1, double>& speeds);
  };
}

#endif // LC_WAVE_PROPAGATION_UPDATER_H

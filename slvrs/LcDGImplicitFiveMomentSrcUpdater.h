/**
 * @file	LcDGImplicitFiveMomentSrcUpdater.h
 *
 * @brief	Implicit updater for 5-moment source terms using a discontinuous Galerkin method
 */

#ifndef LC_DG_IMPLICIT_FIVE_MOMENT_SRC_UPDATER_H
#define LC_DG_IMPLICIT_FIVE_MOMENT_SRC_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Update multi-fluid 5-moment sources using an implicit method. This
 * removes the time-step restriction from plasma- and
 * cyclotron-frequency.
 */
  template <unsigned NDIM>
  class DGImplicitFiveMomentSrcUpdater : public Lucee::UpdaterIfc
  {
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

    private:
/** Pointer to configuration-space basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *basis;
/** Number of fluids */
      unsigned nFluids;
/** Charges of each fluid */
      std::vector<double> charge;
/** Mass of each fluid */
      std::vector<double> mass;
/** Permittivity of free space */
      double epsilon0;
/** Charge-mass ratio for each fluid */
      std::vector<double> qbym;
/** Charge-mass ratio squared for each fluid */
      std::vector<double> qbym2;
/** Flag to indicate if to update pressure equation */
      bool hasPressure;
/** Propagation speed factor for electric field error potential. This
 * is dimensionless and the actual speed is chi_e*c. */
      double chi_e;
/** Propagation speed factor for magnetic field error potential. This
 * is dimensionless and the actual speed is chi_m*c. */
      double chi_m;
/** Electric field error potential damping factor */
      double damp_e;
/** Magnetic field error potential damping factor */
      double damp_m;


  };
}

#endif // LC_IMPLICIT_FIVE_MOMENT_SRC_UPDATER_H

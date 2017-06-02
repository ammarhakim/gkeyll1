/**
 * @file	LcImplicitTenMomentSrcUpdater.h
 *
 * @brief	Implicit updater for 10-moment source terms
 */

#ifndef LC_IMPLICIT_TEN_MOMENT_SRC_UPDATER_H
#define LC_IMPLICIT_TEN_MOMENT_SRC_UPDATER_H

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
 * Update multi-fluid 10-moment sources using an implicit method. This
 * removes the time-step restriction from plasma- and
 * cyclotron-frequency.
 */
  template <unsigned NDIM>
  class ImplicitTenMomentSrcUpdater : public Lucee::UpdaterIfc
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
/** Type of linear solver to use */
      int linSolType;
/** Flag to indicate if static magnetic field is present */
      bool hasStatic;
/** Flag to indicate if sigma is present */
      bool hasSigma; 
/** Fllag to indicate if negative pressure components are restored to previous values */
      bool resetPr;
/** Direction of gravity */
      unsigned grvDir;
/** Gravitational acceleration */
      double gravity;

/**
 * Compute index for fluid current component.
 *
 * @param n Fluid number
 * @param c Component number (0,1,2)
 * @return index for fluid current.
 */
      inline
      int fidx(int n, int c)
      {
        return 3*n+c;
      }

/**
 * Compute index for electric field.
 *
 * @param c Component number (0,1,2)
 * @return index for electric field
 */
      inline
      int eidx(int c)
      {
        return 3*nFluids+c;
      }
  };
}

#endif // LC_IMPLICIT_TEN_MOMENT_SRC_UPDATER_H

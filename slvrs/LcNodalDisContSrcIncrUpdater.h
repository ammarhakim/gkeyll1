/**
 * @file	LcNodalDisContSrcIncrUpdater.h
 *
 * @brief	Updater to compute increment from source terms for use in DG schemes
 */

#ifndef LC_NODAL_DIS_CONT_SRC_INCR_UPDATER_H
#define LC_NODAL_DIS_CONT_SRC_INCR_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcPointSourceIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to integrate ODEs on a grid.
 */
  template <unsigned NDIM>
  class NodalDisContSrcIncrUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new integrator.
 */
      NodalDisContSrcIncrUpdater();

/** Dtor */
      ~NodalDisContSrcIncrUpdater();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

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
/** Point sources to use as RHS in integrator */
      std::vector<Lucee::PointSourceIfc*> rhs;
/** Vector for use in computing sources */
      std::vector<double> ts;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;

/**
 * Compute sources, summing up contributions from each RHS term.
 *
 * @param tm Time at which source is requested.
 * @param xc Coordinates at which source is needed
 * @param inp Inputs for which source is needed
 * @param src On output, sources.
 */
      void calcSource(double tm, const double xc[3], const double *inp, std::vector<double>& src);
  };
}

#endif // LC_NODAL_DIS_CONT_SRC_INCR_UPDATER_H

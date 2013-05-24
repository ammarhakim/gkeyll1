/**
 * @file	LcNodalHyperDiffusionUpdater.h
 *
 * @brief	Updater to evaluate (hyper)diffusion operators using nodal DG
 */

#ifndef LC_NODAL_HYPER_DIFFUSION_UPDATER_H
#define LC_NODAL_HYPER_DIFFUSION_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcStructGridField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to evaluate diffusion (and hyperdiffusion) operator using
 * nodal DG scheme. The stencil used in the update is provided by an
 * instance of the NodalFiniteElementIfc class.
 */
  template <unsigned NDIM>
  class NodalHyperDiffusionUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      NodalHyperDiffusionUpdater();

/**
 * Destroy updater.
 */
      ~NodalHyperDiffusionUpdater();

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
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;
  };
}

#endif // LC_NODAL_HYPER_DIFFUSION_UPDATER_H

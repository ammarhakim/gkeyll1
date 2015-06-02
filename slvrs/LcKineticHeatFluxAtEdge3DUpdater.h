/**
 * @file	LcKineticHeatFluxAtEdge3DUpdater.h
 *
 * @brief	Updater to compute heat flux at right-most edge in domain
 * Used for kinetic SOL problem with 1D2V
 */

#ifndef LC_KINETIC_HEAT_FLUX_AT_EDGE_3D_UPDATER_H
#define LC_KINETIC_HEAT_FLUX_AT_EDGE_3D_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDynVector.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
  class KineticHeatFluxAtEdge3DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      KineticHeatFluxAtEdge3DUpdater();

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
      Lucee::NodalFiniteElementIfc<1> *nodalBasis;
/** Mass of ions in system */
      double ionMass;
/** Mass of electrons in the system */
      double elcMass;
/** Flag specifying whether to compute sheath transmission coefficients also */
      bool computeSheathCoefficient;
  };
}

#endif // LC_KINETIC_HEAT_FLUX_AT_EDGE_3D_UPDATER_H

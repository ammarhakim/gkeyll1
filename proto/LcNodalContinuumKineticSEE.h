/**
 * @file	LcNodalContinuumKineticSEE.h
 *
 * @brief	Updater for calculation Secondary Emitted Electron (SEE) BC
 */

#ifndef LC_NODAL_CONTINUUM_KINETIC_SEE_H
#define LC_NODAL_CONTINUUM_KINETIC_SEE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// math include
#define _USE_MATH_DEFINES
#include <math.h>

namespace Lucee
{
  template <unsigned CDIM, unsigned VDIM>
  class NodalContinuumKineticSEE : public Lucee::UpdaterIfc
  {
    // Number of components for coordinate arrays etc.
    static const unsigned PNC = Lucee::CDIM<CDIM+VDIM>::N;
    static const unsigned CNC = Lucee::CDIM<CDIM>::N;
    public:
      /** Class id: this is used by registration system */
      static const char *id;

      /** Constructor */
      NodalContinuumKineticSEE();
      /** Destructor */
      virtual ~NodalContinuumKineticSEE();
     
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
      // Pointer to phase space basis functions 
      Lucee::NodalFiniteElementIfc<CDIM+VDIM> *phaseBasis;
      // Pointer to configuration space basis functions
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;

      unsigned dir, edge;
  };
}

#endif // LC_NODAL_CONTINUUM_KINETIC_SEE_H

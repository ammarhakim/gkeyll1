/**
 * @file	LcSheathParticleSource1x1v.h
 *
 * @brief	Compute particle source to compensate for loss at boundaries
 */

#ifndef LC_SHEATH_PARTICLE_SOURCE_1X_1V_H
#define LC_SHEATH_PARTICLE_SOURCE_1X_1V_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Updater to compute particle source to compensate for loss of
 * particles at boundaries.
 */
  class SheathParticleSource1x1v : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SheathParticleSource1x1v();

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
      Lucee::NodalFiniteElementIfc<2> *nodalBasis;
/** Thermal velocity of source particles */
      double vThermal;
/** Drift speed of source particles */
      double uDrift;
  };
}

#endif // LC_SHEATH_PARTICLE_SOURCE_1X_1V_H

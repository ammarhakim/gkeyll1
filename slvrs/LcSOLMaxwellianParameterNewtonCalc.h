/**
 * @file  LcSOLMaxwellianParameterNewtonCalc.h
 *
 * @brief	Extremely basic updater that determines next parameters to use as inputs for SOLMaxwellianAtNodeCalc
 */

#ifndef LC_SOL_MAXWELLIAN_PARAMETER_NEWTON_CALC
#define LC_SOL_MAXWELLIAN_PARAMETER_NEWTON_CALC

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includres
#include <Eigen/QR>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class SOLMaxwellianParameterNewtonCalc : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOLMaxwellianParameterNewtonCalc();

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
/** Pointer to configuration space basis functions */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Relative size of forward difference step for u */
      double epsilonU;
/** Relative size of forward difference step for T */
      double epsilonT;
  };
}

#endif // LC_SOL_MAXWELLIAN_PARAMETER_NEWTON_CALC

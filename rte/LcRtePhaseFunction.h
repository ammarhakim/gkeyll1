/**
 * @file	LcRtePhaseFunction.h
 *
 * @brief	Phase function for use in RTE solver.
 */

#ifndef LC_RTE_PHASE_FUNCTION_H
#define LC_RTE_PHASE_FUNCTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcLuaTable.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Base class for phase functions. Children classes should provide a
 * single method that computes and returns the expansion coefficients
 * of the normalized phase function, with normalization \beta_0 = 1.
 */
  class RtePhaseFunction : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Delete the phase-function object */
      virtual ~RtePhaseFunction();

/**
 * Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl) = 0;

/**
 * Return the first L+1 expansion coefficients of the phase
 * function. The first coefficient should be 1.0.
 *
 * @param L number of coefficients required.
 * @return Vector of coefficients. Should be of size L+1.
 */
      virtual Lucee::Vector<double> getExpCoeffs(unsigned L) = 0;

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * This Lua callable method prints out first L+1 phase function
 * coefficients. The lua function should pass a single integer L and
 * has the signature:
 *
 *   phaseFunction:print(L)
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaPrintExpCoeffs(lua_State *L);
  };
}

#endif // LC_RTE_PHASE_FUNCTION_H

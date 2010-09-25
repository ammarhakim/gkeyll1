/**
 * @file	LcRtePhaseFunction.h
 *
 * @brief	Phase function for use in RTE solver.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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
  };
}

#endif // LC_RTE_PHASE_FUNCTION_H

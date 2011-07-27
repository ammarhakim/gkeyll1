/**
 * @file	LcHgPhaseFunction.h
 *
 * @brief	Henyey-Greenstien phase function.
 */

#ifndef LC_HG_PHASE_FUNCTION_H
#define LC_HG_PHASE_FUNCTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaTable.h>
#include <LcRtePhaseFunction.h>

namespace Lucee
{
/**
 * One-term Henyey-Greenstien phase function.
 */
  class HgPhaseFunction : public Lucee::RtePhaseFunction
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Return the first L+1 expansion coefficients of the phase
 * function. The first coefficient should be 1.0.
 *
 * @param L number of coefficients required.
 * @return Vector of coefficients.
 */
      Lucee::Vector<double> getExpCoeffs(unsigned L);

    private:
/** Coefficient of HG phase function */
      double g;
  };
}

#endif // LC_HG_PHASE_FUNCTION_H

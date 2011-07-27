/**
 * @file	LcPlCoeffsPhaseFunction.h
 *
 * @brief	Henyey-Greenstien phase function.
 */

#ifndef LC_PL_COEFFS_PHASE_FUNCTION_H
#define LC_PL_COEFFS_PHASE_FUNCTION_H

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
 * Phase function from input-file specified coefficients.
 */
  class PlCoeffsPhaseFunction : public Lucee::RtePhaseFunction
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
/** Coefficients */
      std::vector<double> coeffs;
  };
}

#endif // LC_PL_COEFFS_PHASE_FUNCTION_H

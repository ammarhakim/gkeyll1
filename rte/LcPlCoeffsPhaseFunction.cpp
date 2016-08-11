/**
 * @file	LcPlCoeffsPhaseFunction.cpp
 *
 * @brief	Henyey-Greenstien phase function.
 */

// lucee includes
#include <LcPlCoeffsPhaseFunction.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *PlCoeffsPhaseFunction::id = "PlCoeffs";

  void
  PlCoeffsPhaseFunction::readInput(Lucee::LuaTable& tbl)
  {
// simply fetch coefficients from table
    coeffs = tbl.getNumVec("coeffs");
// check normalization
    if (coeffs[0] != 1.0)
      throw Lucee::Except(
        "Lucee::PlCoeffsPhaseFunction: Phase function must have 1.0 as first coefficient");
  }
  
  Lucee::Vector<double>
  PlCoeffsPhaseFunction::getExpCoeffs(unsigned L)
  {
    Lucee::Vector<double> pl(L+1);
    unsigned n = (L+1) <= coeffs.size() ? L+1 : coeffs.size();
    for (unsigned i=0; i<n; ++i)
      pl[i] = coeffs[i];
// set remaining to 0.0
    for (unsigned i=n; i<=L; ++i)
      pl[i] = 0.0;

    return pl;
  }
}

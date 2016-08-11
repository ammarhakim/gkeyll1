/**
 * @file	LcHgPhaseFunction.cpp
 *
 * @brief	Henyey-Greenstien phase function.
 */

// lucee includes
#include <LcExcept.h>
#include <LcHgPhaseFunction.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *HgPhaseFunction::id = "HG";

  void
  HgPhaseFunction::readInput(Lucee::LuaTable& tbl)
  {
    if (tbl.hasNumber("g"))
      g = tbl.getNumber("g");
    else
      throw Lucee::Except("HgPhaseFunction::readInput: must provide \"g\" parameter");
  }

  Lucee::Vector<double>
  HgPhaseFunction::getExpCoeffs(unsigned L)
  {
    Lucee::Vector<double> coeffs(L+1);

    coeffs[0] = 1.0; // normalization
    double prod = g;
    for (unsigned i=1; i<=L; ++i)
    {
      coeffs[i] = (2*i+1)*prod;
      prod = prod*g;
    }
    return coeffs;
  }
}

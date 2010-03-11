/**
 * @file	LcHgPhaseFunction.cpp
 *
 * @brief	Henyey-Greenstien phase function.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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

    return coeffs;
  }
}

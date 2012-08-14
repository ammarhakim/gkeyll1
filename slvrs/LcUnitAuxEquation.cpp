/**
 * @file        LcUnitAuxEquation.cpp
 *
 * @brief	Auxilary equations with unit flux.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnitAuxEquation.h>

namespace Lucee
{
// set id for creators
  const char *UnitAuxEquation::id = "Auxiliary";

  UnitAuxEquation::UnitAuxEquation()
    : Lucee::HyperEquation(1, 0)
  {
  }

  void
  UnitAuxEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);

// read number of equations in system
    if (tbl.hasNumber("numEquations"))
      this->setNumEqns((unsigned) tbl.getNumber("numEquations"));
    else
      throw Lucee::Except("UnitAuxEquation::readInput: Must specify number of equations 'numEquations'");

    coeffs.resize(6);
    for (unsigned i=0; i<6; ++i)
      coeffs[i] = 1.0;
// check if coefficients are specified
    if (tbl.hasNumVec("coefficients"))
    {
      std::vector<double> nv = tbl.getNumVec("coefficients");
      for (unsigned i=0; i<nv.size(); ++i)
        coeffs[i] = nv[i];
    }
  }

  void
  UnitAuxEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
  }


  double
  UnitAuxEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    return 0.0;
  }
}

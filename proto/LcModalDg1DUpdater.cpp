/**
 * @file	LcModalDg1DUpdater.cpp
 *
 * @brief	Updater to solver 1D hyperbolic equations using modal DG scheme
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMathLib.h>
#include <LcModalDg1DUpdater.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *ModalDg1DUpdater::id = "ModalDg1D";

  ModalDg1DUpdater::ModalDg1DUpdater()
    : UpdaterIfc(), numBasis(1), Pmk(1,1), DPmk(1,1), normCoeff(1), w(1), mu(1)
  {
  }

  void
  ModalDg1DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// equation to solve
    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("ModalDg1DUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }

// get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL if not explicitly specified
    if (tbl.hasNumber("cflm"))
      cflm = tbl.getNumber("cflm"); // maximum CFL number

// allocate space
    Pmk = Matrix<double>(numBasis, numBasis);
    DPmk = Matrix<double>(numBasis, numBasis);
    normCoeff = Lucee::Vector<double>(numBasis);
    w = Lucee::Vector<double>(numBasis);
    mu = Lucee::Vector<double>(numBasis);

// compute weights and ordinates
    Lucee::gauleg(numBasis, -1, 1, mu, w);

// compute Legendre polynomials at ordinates
    for (unsigned m=0; m<numBasis; ++m)
      for (unsigned k=0; k<numBasis; ++k)
        Pmk(m,k) = Lucee::legendrePoly(m, mu[k]);

// compute derivatives of Legendre polynomials at ordinates
    for (unsigned m=0; m<numBasis; ++m)
      for (unsigned k=0; k<numBasis; ++k)
        DPmk(m,k) = Lucee::legendrePolyDeriv(m, mu[k]);

// compute normalization coefficients
    for (unsigned m=0; m<numBasis; ++m)    
      normCoeff[m] = 1/(2.0*m+1);
  }

  void
  ModalDg1DUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

  }

  Lucee::UpdaterStatus
  ModalDg1DUpdater::update(double t)
  {
    return Lucee::UpdaterStatus();
  }

  void
  ModalDg1DUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}

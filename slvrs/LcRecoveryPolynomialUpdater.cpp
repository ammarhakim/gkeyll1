/**
 * @file	LcRecoveryPolynomialUpdater.cpp
 *
 * @brief	Updater to compute drift velocity and thermal velocity squared
 * given velocity moments 0-2 of the distribution function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcRecoveryPolynomialUpdater.h>
#include <LcStructuredGridBase.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *RecoveryPolynomialUpdater::id = "RecoveryPolynomialUpdater";

  RecoveryPolynomialUpdater::RecoveryPolynomialUpdater()
    : UpdaterIfc()
  {
  }

  void 
  RecoveryPolynomialUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except(
        "RecoveryPolynomialUpdater::readInput: Must specify 2D element to use using 'basis2d'");

    // get hold of 1D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except(
        "RecoveryPolynomialUpdater::readInput: Must specify 1D element to use using 'basis1d'");
  }

  void 
  RecoveryPolynomialUpdater::initialize()
  {
  }

  Lucee::UpdaterStatus 
  RecoveryPolynomialUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Distribution function
    const Lucee::Field<2, double>& fIn = this->getInp<Lucee::Field<2, double> >(0);
    // Recovery polynomial eval on bot edges
    Lucee::Field<1, double>& frBotOut = this->getOut<Lucee::Field<1, double> >(0);
    // Recovery polynomial eval on top edges
    Lucee::Field<1, double>& frTopOut = this->getOut<Lucee::Field<1, double> >(1);

    int nlocal2d = nodalBasis2d->getNumNodes();
    int nlocal1d = nodalBasis1d->getNumNodes();

    if (nlocal2d != 4)
      throw Lucee::Except("RecoveryPolynomialUpdater::update: must use degree 1 basis functions!");

    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> fLeftPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fRightPtr = fIn.createConstPtr();
    Lucee::FieldPtr<double> frTopPtr = frTopOut.createPtr();
    Lucee::FieldPtr<double> frBotPtr = frBotOut.createPtr();

    frTopPtr = 0.0;
    frBotPtr = 0.0;

    // Loop over all x-direction cells
    for (int i = localRgn.getLower(0); i < localRgn.getUpper(0); i++)
    {
      // Bottom edge in v-direction
      int edgeIndex = globalRgn.getLower(1);
      fIn.setPtr(fLeftPtr, i, edgeIndex-1);
      fIn.setPtr(fRightPtr, i, edgeIndex);
      frBotOut.setPtr(frBotPtr, i);

      double fr0 = fLeftPtr[0]/12.0 + fLeftPtr[2]*5.0/12.0 + fRightPtr[0]*5.0/12.0
        + fRightPtr[2]/12.0;

      for (int componentIndex = 0; componentIndex < nlocal1d; componentIndex++)
        frBotPtr[componentIndex] = fr0;

      // Top edge in v-direction
      edgeIndex = globalRgn.getUpper(1) + 1;
      fIn.setPtr(fLeftPtr, i, edgeIndex-1);
      fIn.setPtr(fRightPtr, i, edgeIndex);
      frTopOut.setPtr(frTopPtr, i);

      fr0 = fLeftPtr[0]/12.0 + fLeftPtr[2]*5.0/12.0 + fRightPtr[0]*5.0/12.0
        + fRightPtr[2]/12.0;

      for (int componentIndex = 0; componentIndex < nlocal1d; componentIndex++)
        frTopPtr[componentIndex] = fr0;
    }

    return Lucee::UpdaterStatus();
  }

  void
  RecoveryPolynomialUpdater::declareTypes()
  {
    // takes four inputs (f) 
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // returns two outputs (fr_bot(x), fr_top(x))
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}

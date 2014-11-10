/**
 * @file	LcConstGravitySrcUpdater.cpp
 *
 * @brief	Updater for gravitational source
 */

// lucee includes
#include <LcConstGravitySrcUpdater.h>
#include <LcStructGridField.h>

// eigen inlcudes
#include <Eigen/Eigen>

namespace Lucee
{
// set ids for module system
  const char *ConstGravitySrcUpdater::id = "ConstGravitySrc";

  void
  ConstGravitySrcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    dir = (unsigned) tbl.getNumber("dir");
    gravity = tbl.getNumber("gravity");
  }

  void
  ConstGravitySrcUpdater::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ConstGravitySrcUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    double dt = t-this->getCurrTime();
// get fluid
    Lucee::Field<2, double>& fluid = this->getOut<Lucee::Field<2, double> >(0);
    Lucee::FieldPtr<double> fPtr = fluid.createPtr();

    int dirIdx = dir+1; // offset as first index is rho

    int idx[2];
    Lucee::Region<2, int> localRgn = fluid.getRegion();
    Lucee::RowMajorSequencer<2> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fluid.setPtr(fPtr, idx);

      double rho = fPtr[0];
// old contribution to kinetic energy from momentum in 'dir' direction
      double keold = 0.5*fPtr[dirIdx]*fPtr[dirIdx]/rho;
// update momentum
      fPtr[dirIdx] += gravity*rho*dt;
// now update energy to account for updated momentum
      fPtr[4] = fPtr[4] - keold + 0.5*fPtr[dirIdx]*fPtr[dirIdx]/rho;
    }
    
    return Lucee::UpdaterStatus();
  }
  

  void
  ConstGravitySrcUpdater::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}


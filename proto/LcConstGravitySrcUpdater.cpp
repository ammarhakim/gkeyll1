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
  template <> const char *ConstGravitySrcUpdater<1>::id = "ConstGravitySrc1D";
  template <> const char *ConstGravitySrcUpdater<2>::id = "ConstGravitySrc2D";
  template <> const char *ConstGravitySrcUpdater<3>::id = "ConstGravitySrc3D";

  template <unsigned NDIM>
  void
  ConstGravitySrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    dir = (unsigned) tbl.getNumber("dir");
    gravity = tbl.getNumber("gravity");
  }

  template <unsigned NDIM>
  void
  ConstGravitySrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ConstGravitySrcUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    double dt = t-this->getCurrTime();
// get fluid
    Lucee::Field<NDIM, double>& fluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> fPtr = fluid.createPtr();

    int dirIdx = dir+1; // offset as first index is rho

    int idx[NDIM];
    Lucee::Region<NDIM, int> localRgn = fluid.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
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
  

  template <unsigned NDIM>
  void
  ConstGravitySrcUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ConstGravitySrcUpdater<1>;
  template class ConstGravitySrcUpdater<2>;
  template class ConstGravitySrcUpdater<3>;
}


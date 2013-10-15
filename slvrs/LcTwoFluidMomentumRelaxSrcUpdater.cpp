/**
 * @file	LcTwoFluidMomentumRelaxSrcUpdater.cpp
 *
 * @brief	Updater to apply momentum relaxation from inter-species collisions
 */

// gkeyll includes
#include <LcStructuredGridBase.h>
#include <LcTwoFluidMomentumRelaxSrcUpdater.h>

// std includes
#include <cmath>

namespace Lucee
{

// set ids for module system
  template <> const char *TwoFluidMomentumRelaxSrcUpdater<1>::id = "TwoFluidMomentumRelaxation1D";
  template <> const char *TwoFluidMomentumRelaxSrcUpdater<2>::id = "TwoFluidMomentumRelaxation2D";
  template <> const char *TwoFluidMomentumRelaxSrcUpdater<3>::id = "TwoFluidMomentumRelaxation3D";

  template <unsigned NDIM>
  void
  TwoFluidMomentumRelaxSrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    elcNu = tbl.getNumber("electronCollisionFrequency");
  }

  template <unsigned NDIM>
  void
  TwoFluidMomentumRelaxSrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  TwoFluidMomentumRelaxSrcUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double dt = t-this->getCurrTime();

    Lucee::Field<NDIM, double>& elcFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& ionFluid = this->getOut<Lucee::Field<NDIM, double> >(1);
    Lucee::FieldPtr<double> elcPtr = elcFluid.createPtr();
    Lucee::FieldPtr<double> ionPtr = ionFluid.createPtr();
    int idx[NDIM];

    Lucee::Region<NDIM, int> localRgn = elcFluid.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      elcFluid.setPtr(ionPtr, idx);
      ionFluid.setPtr(elcPtr, idx);
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  TwoFluidMomentumRelaxSrcUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class TwoFluidMomentumRelaxSrcUpdater<1>;
  template class TwoFluidMomentumRelaxSrcUpdater<2>;
  template class TwoFluidMomentumRelaxSrcUpdater<3>;
}


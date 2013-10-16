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
    elcNu = tbl.getNumber("electronIonCollisionFrequency");
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

    double diffU[3];

    Lucee::Region<NDIM, int> localRgn = elcFluid.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      elcFluid.setPtr(ionPtr, idx);
      ionFluid.setPtr(elcPtr, idx);

      double ionNu = elcPtr[0]/ionPtr[0]*elcNu; // ion-electron collision frequency
      double aNu = 0.5*(ionNu+elcNu);
      double ent = std::exp(-aNu*dt);

// initial velocity difference
      for (unsigned d=0; d<3; ++d)
        diffU[d] = elcPtr[1+d]/elcPtr[0]-ionPtr[1+d]/ionPtr[0];

// update electron momentum
      for (unsigned d=0; d<3; ++d)
        elcPtr[1+d] = elcPtr[1+d]-elcNu/(2*aNu)*elcPtr[0]*diffU[d]*(1-ent);
// update ion momentum
      for (unsigned d=0; d<3; ++d)
        ionPtr[1+d] = ionPtr[1+d]+ionNu/(2*aNu)*ionPtr[0]*diffU[d]*(1-ent);
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


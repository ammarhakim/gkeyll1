/**
 * @file	LcFiveMomentNumDensityRelax.cpp
 *
 * @brief	Updater to apply momentum relaxation from inter-species collisions
 */

// gkeyll includes
#include <LcStructuredGridBase.h>
#include <LcFiveMomentNumDensityRelax.h>

// std includes
#include <cmath>

namespace Lucee
{

// set ids for module system
  template <> const char *FiveMomentNumDensityRelax<1>::id = "FiveMomentNumDensityRelax1D";
  template <> const char *FiveMomentNumDensityRelax<2>::id = "FiveMomentNumDensityRelax2D";
  template <> const char *FiveMomentNumDensityRelax<3>::id = "FiveMomentNumDensityRelax3D";

  template <unsigned NDIM>
  void
  FiveMomentNumDensityRelax<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    elcMass = tbl.getNumber("elecronMass");
    ionMass = tbl.getNumber("ionMass");
    gasGamma = tbl.getNumber("gasGamma");
  }

  template <unsigned NDIM>
  void
  FiveMomentNumDensityRelax<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FiveMomentNumDensityRelax<NDIM>::update(double t)
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

      double elcNum0 = elcPtr[0]/elcMass;
      double ionNum0 = ionPtr[0]/ionMass;
      double relaxNum = 0.5*(elcNum0+ionNum0);

// adust electron conserved quantities
      double rho = elcPtr[0];
      double pr = (gasGamma-1)*(elcPtr[4]
        -0.5*(elcPtr[1]*elcPtr[1]+elcPtr[2]*elcPtr[2]+elcPtr[3]*elcPtr[3])/rho);

      elcPtr[0] = relaxNum*elcMass;
      elcPtr[1] = elcPtr[1]/rho*elcPtr[0];
      elcPtr[2] = elcPtr[2]/rho*elcPtr[0];
      elcPtr[3] = elcPtr[3]/rho*elcPtr[0];
      elcPtr[4] = pr/(gasGamma-1) + 0.5*(elcPtr[1]*elcPtr[1]+elcPtr[2]*elcPtr[2]+elcPtr[3]*elcPtr[4])/elcPtr[0];

// adust ion conserved quantities
      rho = ionPtr[0];
      pr = (gasGamma-1)*(ionPtr[4]
        -0.5*(ionPtr[1]*ionPtr[1]+ionPtr[2]*ionPtr[2]+ionPtr[3]*ionPtr[3])/rho);

      ionPtr[0] = relaxNum*ionMass;
      ionPtr[1] = ionPtr[1]/rho*ionPtr[0];
      ionPtr[2] = ionPtr[2]/rho*ionPtr[0];
      ionPtr[3] = ionPtr[3]/rho*ionPtr[0];
      ionPtr[4] = pr/(gasGamma-1) + 0.5*(ionPtr[1]*ionPtr[1]+ionPtr[2]*ionPtr[2]+ionPtr[3]*ionPtr[4])/ionPtr[0];
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  FiveMomentNumDensityRelax<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class FiveMomentNumDensityRelax<1>;
  template class FiveMomentNumDensityRelax<2>;
  template class FiveMomentNumDensityRelax<3>;
}


/**
 * @file	LcConvertTenMomentPtoTUpdater.cpp
 *
 * @brief       Convert ten moment Pressure to temperature
 */

// gkeyll includes
#include <LcConvertTenMomentPToTUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{

// set ids for module system
  template <> const char *ConvertTToPUpdater<1>::id = "ConvertTToP1D";
  template <> const char *ConvertTToPUpdater<2>::id = "ConvertTToP2D";
  template <> const char *ConvertTToPUpdater<3>::id = "ConvertTToP3D";

  template <unsigned NDIM>
  void
  ConvertTToPUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    std::string nf = tbl.getString("mode");
    if (nf == "ToPrimitive"){
      mode = TO_PRIMITIVE;
    } else if (nf == "ToConservative") {
      mode = TO_CONSERVATIVE;
    } else {
      Lucee::Except lce("ConvertTToPUpdater::readInput: 'mode' ");
      lce << nf << " not recognized!" << std::endl;
      throw lce;
    }
    T0 = 0.0;
  }

  template <unsigned NDIM>
  void
  ConvertTToPUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ConvertTToPUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Field<NDIM, double>& tmFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    //    Lucee::Field<NDIM, double>& tMeanFluid = this->getOut<Lucee::Field<NDIM, double> >(1);
    Lucee::FieldPtr<double> ptr = tmFluid.createPtr();
    //    Lucee::FieldPtr<double> meanPtr = tMeanFluid.createPtr();
    int idx[NDIM];

    Lucee::Region<NDIM, int> localRgn = tmFluid.getExtRegion();
    // necessary to convert ghost cells for diffusion
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      tmFluid.setPtr(ptr, idx);
      //      tMeanFluid.setPtr(meanPtr, idx);
      if (mode == TO_PRIMITIVE){
        double r = ptr[0];
        double u = ptr[1]/r;
        double v = ptr[2]/r;
        double w = ptr[3]/r;
        ptr[4] = ptr[4]/r-u*u;
        ptr[5] = ptr[5]/r-u*v;
        ptr[6] = ptr[6]/r-u*w;
        ptr[7] = ptr[7]/r-v*v;
        ptr[8] = ptr[8]/r-v*w;
        ptr[9] = ptr[9]/r-w*w;
        //        T0 = (ptr[4] + ptr[7] + ptr[9])/3.0*0.0;
        //        ptr[4] = ptr[4] - T0;
        //        ptr[7] = ptr[7] - T0;
        //        ptr[9] = ptr[9] - T0;
        //meanPtr[0] = T0;
      } else {
        double r = ptr[0];
        double u = ptr[1]/r;
        double v = ptr[2]/r;
        double w = ptr[3]/r;
        //        T0 = meanPtr[0];
        ptr[4] = (ptr[4])*r+r*u*u;
        ptr[5] = ptr[5]*r+r*u*v;
        ptr[6] = ptr[6]*r+r*u*w;
        ptr[7] = (ptr[7])*r+r*v*v;
        ptr[8] = ptr[8]*r+r*v*w;
        ptr[9] = (ptr[9])*r+r*w*w;
      }
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ConvertTToPUpdater<NDIM>::declareTypes()
  {
    //this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
        this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ConvertTToPUpdater<1>;
  template class ConvertTToPUpdater<2>;
  template class ConvertTToPUpdater<3>;
}


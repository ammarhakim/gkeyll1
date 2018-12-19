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
  template <> const char *ConvertPToTUpdater<1>::id = "ConvertPToT1D";
  template <> const char *ConvertPToTUpdater<2>::id = "ConvertPToT2D";
  template <> const char *ConvertPToTUpdater<3>::id = "ConvertPToT3D";

  template <unsigned NDIM>
  void
  ConvertPToTUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    std::string nf = tbl.getString("mode");
    if (nf == "ToPrimitive"){
      mode = TO_PRIMITIVE;
    } else if (nf == "ToConservative") {
      mode = TO_CONSERVATIVE;
    } else {
      Lucee::Except lce("ConvertPToTUpdater::readInput: 'mode' ");
      lce << nf << " not recognized!" << std::endl;
      throw lce;
    }

    outType = OUT_TEMPERATURE;
    if (tbl.hasString("outType")) { 
      std::string outQty = tbl.getString("outType");
      if (outQty == "OutPressure"){
        outType = OUT_PRESSURE;
      } else if (outQty == "OutTemperature") {
        outType = OUT_TEMPERATURE;
      } else {
        Lucee::Except lce("ConvertPToTUpdater::readInput: 'outType' ");
        lce << outQty << " not recognized!" << std::endl;
        throw lce;
      }
    }
  }

  template <unsigned NDIM>
  void
  ConvertPToTUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ConvertPToTUpdater<NDIM>::update(double t)
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
        // this converts pressure to temperature
        if (outType == OUT_TEMPERATURE) {
          ptr[4] = ptr[4]/r-u*u;
          ptr[5] = ptr[5]/r-u*v;
          ptr[6] = ptr[6]/r-u*w;
          ptr[7] = ptr[7]/r-v*v;
          ptr[8] = ptr[8]/r-v*w;
          ptr[9] = ptr[9]/r-w*w;
        } else {
          // convert to pressure
          ptr[4] = ptr[4]-u*u*r;
          ptr[5] = ptr[5]-u*v*r;
          ptr[6] = ptr[6]-u*w*r;
          ptr[7] = ptr[7]-v*v*r;
          ptr[8] = ptr[8]-v*w*r;
          ptr[9] = ptr[9]-w*w*r;       
        }
      } else {
        double r = ptr[0];
        double u = ptr[1]/r;
        double v = ptr[2]/r;
        double w = ptr[3]/r;
        //        T0 = meanPtr[0];
        if (outType == OUT_TEMPERATURE) {
          // given temperature, convert to P + rho u u
          ptr[4] = ptr[4]*r+r*u*u;
          ptr[5] = ptr[5]*r+r*u*v;
          ptr[6] = ptr[6]*r+r*u*w;
          ptr[7] = ptr[7]*r+r*v*v;
          ptr[8] = ptr[8]*r+r*v*w;
          ptr[9] = ptr[9]*r+r*w*w;
        } else {
          // give P, convert to P + rho u u 
          ptr[4] = ptr[4] + r*u*u;
          ptr[5] = ptr[5] + r*u*v;
          ptr[6] = ptr[6] + r*u*w;
          ptr[7] = ptr[7] + r*v*v;
          ptr[8] = ptr[8] + r*v*w;
          ptr[9] = ptr[9] + r*w*w;
        }
      }
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ConvertPToTUpdater<NDIM>::declareTypes()
  {
    //this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
        this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ConvertPToTUpdater<1>;
  template class ConvertPToTUpdater<2>;
  template class ConvertPToTUpdater<3>;
}


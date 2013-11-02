/**
 * @file	LcTenMomLocalAnisoHeatFluxUpdater.cpp
 *
 * @brief	Implicit updater for 10-moment collisional source terms
 */

// gkeyll includes
#include <LcTenMomLocalAnisoHeatFluxUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{

// set ids for module system
  template <> const char *TenMomLocalAnisoHeatFluxUpdater<1>::id = "TenMomentLocalAnisoCollisionlessHeatFlux1D";
  template <> const char *TenMomLocalAnisoHeatFluxUpdater<2>::id = "TenMomentLocalAnisoCollisionlessHeatFlux2D";
  template <> const char *TenMomLocalAnisoHeatFluxUpdater<3>::id = "TenMomentLocalAnisoCollisionlessHeatFlux3D";

  template <unsigned NDIM>
  void
  TenMomLocalAnisoHeatFluxUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    kA = tbl.getNumber("averageWaveNumber");
  }

  template <unsigned NDIM>
  void
  TenMomLocalAnisoHeatFluxUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  TenMomLocalAnisoHeatFluxUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double dt = t-this->getCurrTime();

    Lucee::Field<NDIM, double>& tmFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> ptr = tmFluid.createPtr();
    int idx[NDIM];

    Lucee::Region<NDIM, int> localRgn = tmFluid.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      tmFluid.setPtr(ptr, idx);

      double r = ptr[0];
      double u = ptr[1]/r;
      double v = ptr[2]/r;
      double w = ptr[3]/r;
      double pxx = ptr[4]-r*u*u;
      double pxy = ptr[5]-r*u*v;
      double pxz = ptr[6]-r*u*w;
      double pyy = ptr[7]-r*v*v;
      double pyz = ptr[8]-r*v*w;
      double pzz = ptr[9]-r*w*w;

      double p = (pxx+pyy+pzz)/3.0;
      double vt = std::sqrt(p/r);
      double edt = std::exp(-vt*kA*dt);

// compute updated pressure tensor component
      ptr[4] = (pxx-p)*edt+p + r*u*u;
      ptr[5] = pxy*edt + r*u*v;
      ptr[6] = pxz*edt + r*u*w;
      ptr[7] = (pyy-p)*edt+p + r*v*v;
      ptr[8] = pyz*edt + r*v*w;
      ptr[9] = (pzz-p)*edt+p + r*w*w;
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  TenMomLocalAnisoHeatFluxUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class TenMomLocalAnisoHeatFluxUpdater<1>;
  template class TenMomLocalAnisoHeatFluxUpdater<2>;
  template class TenMomLocalAnisoHeatFluxUpdater<3>;
}


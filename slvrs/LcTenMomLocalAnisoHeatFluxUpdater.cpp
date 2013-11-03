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
// indices into symmetric 3x3 matrix
  static const unsigned T11 = 0;
  static const unsigned T12 = 1;
  static const unsigned T13 = 2;
  static const unsigned T22 = 3;
  static const unsigned T23 = 4;
  static const unsigned T33 = 5;
// indices for magentic field
  static const unsigned BX = 3;
  static const unsigned BY = 4;
  static const unsigned BZ = 5;

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
    B0 = tbl.getNumber("refMagneticField");
// this may need more clever selection
    eps = 1e-14;
    if (B0>0)
      eps = B0*1e-6;
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

    const Lucee::Field<NDIM, double>& emField = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& tmFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> ptr = tmFluid.createPtr();
    Lucee::ConstFieldPtr<double> emPtr = tmFluid.createConstPtr();
    int idx[NDIM];

    double N[6], T[6], cgl[6], p[6], u[3];
    double e3 = eps/3.0;

    Lucee::Region<NDIM, int> localRgn = tmFluid.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      tmFluid.setPtr(ptr, idx);
      emField.setPtr(emPtr, idx);

      double bNorm = emPtr[BX]*emPtr[BX]+emPtr[BY]*emPtr[BY]+emPtr[BZ]*emPtr[BZ];
// construct needed tensors: N projects pressure to parallel
// direction, T to perpendicular direction. The eps terms ensure that
// in the absence of a magnetic field the relaxation is isotropic.
      N[T11] = (emPtr[BX]*emPtr[BX] + e3)/(bNorm+eps);
      N[T12] = emPtr[BX]*emPtr[BY]/(bNorm+eps);
      N[T13] = emPtr[BX]*emPtr[BZ]/(bNorm+eps);
      N[T22] = (emPtr[BY]*emPtr[BY] + e3)/(bNorm+eps);
      N[T23] = emPtr[BY]*emPtr[BZ]/(bNorm+eps);
      N[T33] = (emPtr[BZ]*emPtr[BZ] + e3)/(bNorm+eps);

      T[T11] = 1-N[T11];
      T[T12] = -N[T12];
      T[T13] = -N[T13];
      T[T22] = 1-N[T22];
      T[T23] = -N[T23];
      T[T33] = 1-N[T33];

      double r = ptr[0];
      u[0] = ptr[1]/r;
      u[1] = ptr[2]/r;
      u[2] = ptr[3]/r;
      p[T11] = ptr[4]-r*u[0]*u[0];
      p[T12] = ptr[5]-r*u[0]*u[1];
      p[T13] = ptr[6]-r*u[0]*u[2];
      p[T22] = ptr[7]-r*u[1]*u[1];
      p[T23] = ptr[8]-r*u[1]*u[2];
      p[T33] = ptr[9]-r*u[2]*u[2];

      double pr = (p[T11]+p[T22]+p[T33])/3;
      double vt = std::sqrt(pr/r);
      double edt = std::exp(-vt*kA*dt);

// compute p_par and p_perp
      double p_par = p[T11]*N[T11]+p[T22]*N[T22]+p[T33]*N[T33]
        +2*(p[T12]*N[T12]+p[T13]*N[T13]+p[T23]*N[T23]);
      double p_per = 0.5*(3*pr-p_par);

// construct CGL part of pressure tensor (this remains invariant under
// the relaxation process)
      for (unsigned i=0; i<6; ++i)
        cgl[i] = p_par*N[i]+p_per*T[i];

// compute updated pressure tensor components
      ptr[4] = (p[T11]-cgl[T11])*edt + cgl[T11] + r*u[0]*u[0];
      ptr[5] = (p[T12]-cgl[T12])*edt + cgl[T12] + r*u[0]*u[1];
      ptr[6] = (p[T13]-cgl[T13])*edt + cgl[T13] + r*u[0]*u[2];
      ptr[7] = (p[T22]-cgl[T22])*edt + cgl[T22] + r*u[1]*u[1];
      ptr[8] = (p[T23]-cgl[T23])*edt + cgl[T23] + r*u[1]*u[2];
      ptr[9] = (p[T33]-cgl[T33])*edt + cgl[T33] + r*u[2]*u[2];
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  TenMomLocalAnisoHeatFluxUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class TenMomLocalAnisoHeatFluxUpdater<1>;
  template class TenMomLocalAnisoHeatFluxUpdater<2>;
  template class TenMomLocalAnisoHeatFluxUpdater<3>;
}


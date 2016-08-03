/**
 * @file	LcImplicitTwentyMomentHeatFluxLimitUpdater.cpp
 *
 * @brief	Implicit updater for 20-moment heat flux limiting
 */

// gkeyll includes
#include <LcImplicitTwentyMomentHeatFluxLimitUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{
  static const unsigned X = 0;
  static const unsigned Y = 1;
  static const unsigned Z = 2;

  static const unsigned RHO = 0;
  static const unsigned RHOUX = 1;
  static const unsigned RHOUY = 2;
  static const unsigned RHOUZ = 3;

  static const unsigned P11 = 4;
  static const unsigned P12 = 5;
  static const unsigned P13 = 6;
  static const unsigned P22 = 7;
  static const unsigned P23 = 8;
  static const unsigned P33 = 9;

  static const unsigned Q111 = 10;
  static const unsigned Q112 = 11;
  static const unsigned Q113 = 12;
  static const unsigned Q122 = 13;
  static const unsigned Q123 = 14;
  static const unsigned Q133 = 15;
  static const unsigned Q222 = 16;
  static const unsigned Q223 = 17;
  static const unsigned Q233 = 18;
  static const unsigned Q333 = 19;

// set ids for module system
  template <> const char *ImplicitTwentyMomentHeatFluxLimitUpdater<1>::id = "ImplicitTwentyMomentHeatFluxLimit1D";
  template <> const char *ImplicitTwentyMomentHeatFluxLimitUpdater<2>::id = "ImplicitTwentyMomentHeatFluxLimit2D";
  template <> const char *ImplicitTwentyMomentHeatFluxLimitUpdater<3>::id = "ImplicitTwentyMomentHeatFluxLimit3D";

  template <unsigned NDIM>
  void
  ImplicitTwentyMomentHeatFluxLimitUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    if (tbl.hasNumber("frac")){
      frac = tbl.getNumber("frac");
    } else {
      frac = 1.0;
    }
    if (tbl.hasNumber("dampingFrac")) {
      dampingFrac = tbl.getNumber("dampingFrac");
    } else { 
      dampingFrac = 0.5;
    }
    if (tbl.hasBool("conservative")){
      conservative = tbl.getBool("conservative");
    } else {
      conservative = true;
    }

  }

  template <unsigned NDIM>
  void
  ImplicitTwentyMomentHeatFluxLimitUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ImplicitTwentyMomentHeatFluxLimitUpdater<NDIM>::update(double t)
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
      double qxxx, qxxy, qxxz, qxyy, qxyz, qxzz, qyyy, qyyz, qyzz, qzzz;
      if (conservative){
        qxxx = ptr[Q111] - 3*pxx*u - r*u*u*u;
        qxxy = ptr[Q112] - 2*pxy*u - pxx*v - r*u*u*v;
        qxxz = ptr[Q113] - 2*pxz*u - pxx*w - r*u*u*w;
        qxyy = ptr[Q122] - pyy*u - 2*pxy*v - r*u*v*v;
        qxyz = ptr[Q123] - pyz*u - pxz*v - pxy*w - r*u*v*w;
        qxzz = ptr[Q133] - pzz*u - 2*pxz*w - r*u*w*w;
        qyyy = ptr[Q222] - 3*pyy*v - r*v*v*v;
        qyyz = ptr[Q223] - 2*pyz*v - pyy*w - r*v*v*w;
        qyzz = ptr[Q233] - pzz*v - 2*pyz*w - r*v*w*w;
        qzzz = ptr[Q333] - 3*pzz*w - r*w*w*w;
      } else {
        qxxx = ptr[Q111];
        qxxy = ptr[Q112];
        qxxz = ptr[Q113];
        qxyy = ptr[Q122];
        qxyz = ptr[Q123];
        qxzz = ptr[Q133];
        qyyy = ptr[Q222];
        qyyz = ptr[Q223];
        qyzz = ptr[Q233];
        qzzz = ptr[Q333];
      }
      // here we throw in a heat flux limiter. 
      // in this form of the operator we will ignore the k vt term
      double reduction = 1.0; 

      // this is the physical realisability limit
      // Here we have R_iijj > Q_kii P^-1 kl Qljj + P_ii P_jj/rho
      //      double r = (3*pxx*pxx + 4*pxy*pxy + 4*pxz*pxz + 3*pyy*pyy + 4*pyz*pyz + 2*pyy*pzz + 3*pzz*pzz + 2*pxx*(pyy + pzz))/p0;
      double rhat = (2*(pxx*pxx + 2*pxy*pxy + 2*pxz*pxz + pyy*pyy + 2*pyz*pyz + pzz*pzz))/r; // rhat = R_iijj - p_ii p_jj/rho
      double q1 = qxxx + qxyy + qxzz;
      double q2 = qyyy + qxxy + qyzz;
      double q3 = qzzz + qxxz + qyyz;
      double det = -pxz*pxz*pyy + 2*pxy*pxz*pyz -pxx*pyz*pyz - pxy*pxy*pzz + pxx*pyy*pzz;
      double PInv[6]; // need to calculate inverse pressure tensor
      PInv[0] = (-pyz*pyz + pyy*pzz)/det;
      PInv[1] = (pxz*pyz - pxy*pzz)/det;
      PInv[2] = (-pxz*pyy + pxy*pyz)/det;
      PInv[3] = (-pxz*pxz + pxx*pzz)/det;
      PInv[4] = (pxy*pxz - pxx*pyz)/det;
      PInv[5] = (-pxy*pxy + pxx*pyy)/det; 
      double qpq = q1*q1*PInv[0] + 2*q1*q2*PInv[1] + 2*q1*q3*PInv[2] + q2*q2*PInv[3] +  2*q2*q3*PInv[4] + q3*q3*PInv[5];
      
      // want qpq < rhat.
      if (qpq > rhat) { 
        reduction = std::sqrt(fabs(rhat/qpq)); 
        std::cout<<" applied reduction of " << reduction <<"\n";
        //        std::cout<<"rhat = " << rhat <<" and qpq = " << qpq << "\n";
      }

      // the heat flux with the largest deviation is reduced to a base value
     
      // new heat flux = (qiii - sign(x)*frac*dampingFrac*sqrt(pii*pii*pii/r))*edt
      // fraction = new heat flux/qiii
      // apply this
      if (conservative) {
        ptr[Q111] = reduction*qxxx + 3*pxx*u + r*u*u*u;
        ptr[Q222] = reduction*qyyy + 3*pyy*v + r*v*v*v;
        ptr[Q333] = reduction*qzzz + 3*pzz*w + r*w*w*w;
        ptr[Q112] = reduction*qxxy + 2*pxy*u + pxx*v + r*u*u*v;
        ptr[Q113] = reduction*qxxz + 2*pxz*u + pxx*w + r*u*u*w;
        ptr[Q122] = reduction*qxyy + pyy*u + 2*pxy*v + r*u*v*v;
        ptr[Q123] = reduction*qxyz + pyz*u + pxz*v + pxy*w + r*u*v*w;
        ptr[Q133] = reduction*qxzz + pzz*u + 2*pxz*w + r*u*w*w;
        ptr[Q223] = reduction*qyyz + 2*pyz*v + pyy*w + r*v*v*w;
        ptr[Q233] = reduction*qyzz + pzz*v + 2*pyz*w + r*v*w*w;
      } else {
        ptr[Q111] = reduction*qxxx;
        ptr[Q222] = reduction*qyyy;
        ptr[Q333] = reduction*qzzz;
        ptr[Q112] = reduction*qxxy;
        ptr[Q113] = reduction*qxxz;
        ptr[Q122] = reduction*qxyy;
        ptr[Q123] = reduction*qxyz;
        ptr[Q133] = reduction*qxzz;
        ptr[Q223] = reduction*qyyz;
        ptr[Q233] = reduction*qyzz;
      }
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ImplicitTwentyMomentHeatFluxLimitUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ImplicitTwentyMomentHeatFluxLimitUpdater<1>;
  template class ImplicitTwentyMomentHeatFluxLimitUpdater<2>;
  template class ImplicitTwentyMomentHeatFluxLimitUpdater<3>;
}


/**
 * @file	LcImplicitTwentyMomentCollisionUpdater.cpp
 *
 * @brief	Implicit updater for 20-moment collisional source terms
 */

// gkeyll includes
#include <LcImplicitTwentyMomentCollisionUpdater.h>
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
  template <> const char *ImplicitTwentyMomentCollisionUpdater<1>::id = "ImplicitTwentyMomentCollisionSrc1D";
  template <> const char *ImplicitTwentyMomentCollisionUpdater<2>::id = "ImplicitTwentyMomentCollisionSrc2D";
  template <> const char *ImplicitTwentyMomentCollisionUpdater<3>::id = "ImplicitTwentyMomentCollisionSrc3D";

  template <unsigned NDIM>
  void
  ImplicitTwentyMomentCollisionUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    nu = tbl.getNumber("collisionFrequency");
    if (tbl.hasBool("qLim")){
        lim = tbl.getBool("qLim");
    } else {
        lim = false;
    }
    if (tbl.hasNumber("frac")){
        frac = tbl.getNumber("frac");
    } else {
        frac = 1.0;
    }

  }

  template <unsigned NDIM>
  void
  ImplicitTwentyMomentCollisionUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ImplicitTwentyMomentCollisionUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double dt = t-this->getCurrTime();

    double edt = std::exp(-nu*dt);

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

      double qxxx = ptr[Q111] - 3*pxx*u - r*u*u*u;
      double qxxy = ptr[Q112] - 2*pxy*u - pxx*v - r*u*u*v;
      double qxxz = ptr[Q113] - 2*pxz*u - pxx*w - r*u*u*w;
      double qxyy = ptr[Q122] - pyy*u - 2*pxy*v - r*u*v*v;
      double qxyz = ptr[Q123] - pyz*u - pxz*v - pxy*w - r*u*v*w;
      double qxzz = ptr[Q133] - pzz*u - 2*pxz*w - r*u*w*w;
      double qyyy = ptr[Q222] - 3*pyy*v - r*v*v*v;
      double qyyz = ptr[Q223] - 2*pyz*v - pyy*w - r*v*v*w;
      double qyzz = ptr[Q233] - pzz*v - 2*pyz*w - r*v*w*w;
      double qzzz = ptr[Q333] - 3*pzz*w - r*w*w*w;

// compute updated pressure tensor component
      ptr[4] = (pxx-p)*edt+p + r*u*u;
      ptr[5] = pxy*edt + r*u*v;
      ptr[6] = pxz*edt + r*u*w;
      ptr[7] = (pyy-p)*edt+p + r*v*v;
      ptr[8] = pyz*edt + r*v*w;
      ptr[9] = (pzz-p)*edt+p + r*w*w;

      ptr[Q111] = edt*qxxx + 3*pxx*u + r*u*u*u;
      ptr[Q112] = edt*qxxy + 2*pxy*u + pxx*v + r*u*u*v;
      ptr[Q113] = edt*qxxz + 2*pxz*u + pxx*w + r*u*u*w;
      ptr[Q122] = edt*qxyy + pyy*u + 2*pxy*v + r*u*v*v;
      ptr[Q123] = edt*qxyz + pyz*u + pxz*v + pxy*w + r*u*v*w;
      ptr[Q133] = edt*qxzz + pzz*u + 2*pxz*w + r*u*w*w;
      ptr[Q222] = edt*qyyy + 3*pyy*v + r*v*v*v;
      ptr[Q223] = edt*qyyz + 2*pyz*v + pyy*w + r*v*v*w;
      ptr[Q233] = edt*qyzz + pzz*v + 2*pyz*w + r*v*w*w;
      ptr[Q333] = edt*qzzz + 3*pzz*w + r*w*w*w;

      if (lim){

        double boundx = sqrt(r)*qxxx/sqrt(pxx)/pxx;
        double boundy = sqrt(r)*qyyy/sqrt(pyy)/pyy;
        double boundz = sqrt(r)*qzzz/sqrt(pzz)/pzz;
        double boundmax = std::max(std::max(fabs(boundx),fabs(boundy)),fabs(boundz));

        if (boundmax > frac) { 
          // set offending q = 1/normalised_q*frac
          double reduction = frac/boundmax;
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
          
        }        
      } 
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ImplicitTwentyMomentCollisionUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ImplicitTwentyMomentCollisionUpdater<1>;
  template class ImplicitTwentyMomentCollisionUpdater<2>;
  template class ImplicitTwentyMomentCollisionUpdater<3>;
}


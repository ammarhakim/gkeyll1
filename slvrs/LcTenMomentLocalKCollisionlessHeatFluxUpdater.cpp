/**
 * @file	LcTenMomentLocalKCollisionlessHeatFluxUpdater.cpp
 *
 * @brief	Implicit updater for 10-moment collisional source terms
 */

// gkeyll includes
#include <LcTenMomentLocalKCollisionlessHeatFluxUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{

// set ids for module system
  template <> const char *TenMomentLocalKCollisionlessHeatFluxUpdater<1>::id = "TenMomentLocalKCollisionlessHeatFlux1D";
  template <> const char *TenMomentLocalKCollisionlessHeatFluxUpdater<2>::id = "TenMomentLocalKCollisionlessHeatFlux2D";
  template <> const char *TenMomentLocalKCollisionlessHeatFluxUpdater<3>::id = "TenMomentLocalKCollisionlessHeatFlux3D";

  template <unsigned NDIM>
  void
  TenMomentLocalKCollisionlessHeatFluxUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    kA = -1.;
    if (tbl.hasNumber("averageWaveNumber"))
      kA = tbl.getNumber("averageWaveNumber");

    kmax = -1.;
    if (tbl.hasNumber("maxWaveNumber"))
      kmax = tbl.getNumber("maxWaveNumber");

    kmin = -1.;
    if (tbl.hasNumber("minWaveNumber"))
      kmin = tbl.getNumber("minWaveNumber");
  }

  template <unsigned NDIM>
  void
  TenMomentLocalKCollisionlessHeatFluxUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  TenMomentLocalKCollisionlessHeatFluxUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double dt = t-this->getCurrTime();

    Lucee::Field<NDIM, double>& tmFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> ptr = tmFluid.createPtr();
    Lucee::FieldPtr<double> rPtr = tmFluid.createPtr();
    Lucee::FieldPtr<double> lPtr = tmFluid.createPtr();
    int idx[NDIM];

    Lucee::Region<NDIM, int> localRgn = tmFluid.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    double k_max_used = 0;
    double k_min_used = 1e10;
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      tmFluid.setPtr(ptr, idx);
      grid.setIndex(idx);

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

      // k for current cell, default is kA if kA > 0
      double kThis = kA;
      // compute kA based on local length scale
      if (kA < 0.)
      {
        // k vector
        double k[3];
        // compoonents used to be compute local length scale
        unsigned components[3];
        for (unsigned d = 0; d < 3; ++d)
        {
          k[d] = 0.;
          components[d] = 0; // FIXME using rho to compute gradient for all direction/tensor component
        }

        // compute k[d]   
        for (unsigned d = 0; d < NDIM; ++d)
        {
          unsigned c = components[d];

          // FIXME ghost cell might not have a neighbor cell
          idx[d] = idx[d]+1;
          tmFluid.setPtr(rPtr, idx); // right cell
          grid.setIndex(idx);
          double xcr[3];
          grid.getCentroid(xcr);
        
          idx[d] = idx[d]-1-1;
          tmFluid.setPtr(lPtr, idx); // left cell
          grid.setIndex(idx);
          double xcl[3];
          grid.getCentroid(xcl);
        
          double dx1 = 1./(xcr[d]-xcl[d]);
          idx[d] = idx[d]+1; // restore idx
          
          double k_d = std::fabs(dx1*(rPtr[c] - lPtr[c])/ptr[c]);

          if (kmax > 0. && k_d > kmax)
            k_d = kmax; // TODO: physics based kmax (e.g., 1/d_e)
          if (kmin > 0. && k_d < kmin)
            k_d = kmin; // TODO: physics based kmin (e.g., wc/vt = qB/m/vt)
        
          k[d] = k_d;
        }
        
        // compute k = |k|
        double k2 = 0.;
        for (unsigned d = 0; d < NDIM; ++d)
          k2 += k[d]*k[d];
        kThis = std::sqrt(k2);
      }
      if (kThis > k_max_used) k_max_used = kThis;
      if (kThis < k_min_used) k_min_used = kThis;

      double edt = std::exp(-vt*kThis*dt);

// compute updated pressure tensor component
      ptr[4] = (pxx-p)*edt+p + r*u*u;
      ptr[5] = pxy*edt + r*u*v;
      ptr[6] = pxz*edt + r*u*w;
      ptr[7] = (pyy-p)*edt+p + r*v*v;
      ptr[8] = pyz*edt + r*v*w;
      ptr[9] = (pzz-p)*edt+p + r*w*w;
    }
    std::cout << "k range [" << k_min_used << ", " << k_max_used << "]" << std::endl;
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  TenMomentLocalKCollisionlessHeatFluxUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class TenMomentLocalKCollisionlessHeatFluxUpdater<1>;
  template class TenMomentLocalKCollisionlessHeatFluxUpdater<2>;
  template class TenMomentLocalKCollisionlessHeatFluxUpdater<3>;
}


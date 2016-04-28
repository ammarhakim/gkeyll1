/**
 * @file	LcPeriodicCollisionlessHeatFluxUpdater.cpp
 *
 * @brief	Updater for Hammett-Perkins with Periodic BCs
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMathPhysConstants.h>
#include <LcPeriodicCollisionlessHeatFluxUpdater.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

namespace Lucee
{
  static const unsigned T11 = 4;
  static const unsigned T12 = 5;
  static const unsigned T13 = 6;
  static const unsigned T22 = 7;
  static const unsigned T23 = 8;
  static const unsigned T33 = 9;

/** Class id: this is used by registration system */
  template <> const char *PeriodicCollisionlessHeatFluxUpdater<1>::id = "PeriodicCollisionlessHeatFlux1D";
  template <> const char *PeriodicCollisionlessHeatFluxUpdater<2>::id = "PeriodicCollisionlessHeatFlux2D";
  template <> const char *PeriodicCollisionlessHeatFluxUpdater<3>::id = "PeriodicCollisionlessHeatFlux3D";

  template <unsigned NDIM>
  PeriodicCollisionlessHeatFluxUpdater<NDIM>::PeriodicCollisionlessHeatFluxUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  PeriodicCollisionlessHeatFluxUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  PeriodicCollisionlessHeatFluxUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();

// get hold of grid
    const Lucee::RectCartGrid<NDIM>& grid 
      = this->getGrid<Lucee::RectCartGrid<NDIM> >();
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    int vol = localRgn.getVolume();

// computational space
    Lucee::Region<NDIM, double> compRgn = grid.getComputationalSpace();
    double Lx = compRgn.getShape(0);
    double Ly = compRgn.getShape(1);
// construct arrays for use in FFTW -- no local allocation in the serial version
// using the FFTW many interface so allocate memory for all 6 arrays.
    src_in_out.resize(vol*6);
    sol_in_out.resize(vol*6);

// construct wave-number arrays. Funkyness is needed as DC component
// is stored in location [0].
    int nx = localRgn.getShape(0);
    kx.resize(nx);
    for (unsigned i=0; i<nx/2; ++i)
      kx[i] = 2*Lucee::PI/Lx*i;
    for (unsigned i=nx/2; i>0; --i)
      kx[nx-i] = -2*Lucee::PI/Lx*i;
    
    int ny = localRgn.getShape(1);
    ky.resize(ny);
    for (unsigned i=0; i<ny/2; ++i)
      ky[i] = 2*Lucee::PI/Ly*i;
    for (unsigned i=ny/2; i>0; --i)
      ky[ny-i] = -2*Lucee::PI/Ly*i;

    kabs.resize(nx*ny); 
    for (unsigned i=0; i<nx; ++i){
      for (unsigned j=0; j<ny; ++j){
        kabs[j+ny*i] = std::sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
      }
    }

// no need to reset DC component since it's a multiply by |k|

// these ugly casts are needed so we can pass pointers of correct type
// to FFTW. This works because for *most* C++ compilers
// complex<double> is the same as double [2] which FFTW uses and is
// typedef-ed as fftw_complex. It is possible that this can cause
// problems if some C++ compiler does not support this.
    f_src_in_out = reinterpret_cast<fftw_complex*> (&src_in_out[0]);
    f_sol_in_out = reinterpret_cast<fftw_complex*> (&sol_in_out[0]);

// now create plans for forward and inverse FFTs
// fft_plan_many_dft(int rank, const int *n, int howmany, fftw_complex *in, const int *inembed, int istride, int idist, fftw_complex * out, const int *oembed, int ostride, int odist, int sign, unsigned flags); 
    int n[2] = {nx,ny}; // size of domain -- FIXME: generalise to NDIM
    int howmany = 6; // transform all 6 components
    int istride = 1; // distance between elements in memory. We could use this to interleave the data for better memory access? 
    int ostride = 1; 
    int idist = nx*ny; // distance between first element of array
    int odist = nx*ny; 
    f_plan = fftw_plan_many_dft(
             2, n, howmany, f_src_in_out, NULL, istride, idist, f_src_in_out, NULL, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE); //could replace estimate by patient?

    b_plan = fftw_plan_many_dft(
             2, n, howmany, f_sol_in_out, NULL, istride, idist, f_sol_in_out, NULL, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  PeriodicCollisionlessHeatFluxUpdater<NDIM>::update(double t)
  {
// get hold of grid and find dt
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double dt = t-this->getCurrTime();
// get fields
    Lucee::Field<NDIM, double>& tmFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> ptr = tmFluid.createPtr();
    int idx[NDIM]; 
    double vTarr[6]; // we are actually calculating T/m, not T
// local indexing region
    Lucee::Region<NDIM, int> localRgn = tmFluid.getRegion();
    Lucee::RowMajorSequencer <NDIM> seq(localRgn);
    Lucee::RowMajorIndexer<NDIM> idxr(localRgn);
    int vol = localRgn.getVolume();
    // total isotropic temperature = Sum (Txx+Tyy+Tzz)/3.0
    double intIsovT = 0.0;
    double intvTxx = 0.0; 
    double intvTyy = 0.0;
    double intvTzz = 0.0;
    // do stepping
    while (seq.step()) { 

      seq.fillWithIndex(idx);
      tmFluid.setPtr(ptr, idx); 
      double r = ptr[0];
      double u = ptr[1]/r;
      double v = ptr[2]/r;
      double w = ptr[3]/r;
      vTarr[0] = (ptr[4]-r*u*u)/r;
      vTarr[1] = (ptr[5]-r*u*v)/r;
      vTarr[2] = (ptr[6]-r*u*w)/r;
      vTarr[3] = (ptr[7]-r*v*v)/r;
      vTarr[4] = (ptr[8]-r*v*w)/r;
      vTarr[5] = (ptr[9]-r*w*w)/r;
      // find global isotropic temperature and fill matrix
      intIsovT += (vTarr[0]+vTarr[3]+vTarr[5])/3.0;
      intvTxx += vTarr[0];
      intvTyy += vTarr[3];
      intvTzz += vTarr[5];
      // fill the matrix
      for (unsigned i = 0; i < 6; ++i){
        src_in_out[idxr.getIndex(idx) + i*vol] = vTarr[i];
      }
    }

    double avgvT = intIsovT/vol; // MPI Allreduce in parallel version
    double avgvTxx = intvTxx/vol; 
    double avgvTyy = intvTyy/vol;
    double avgvTzz = intvTzz/vol;
    // subtract isotopric temperature from input data
    seq.reset(); // reset the sequencer
      
    while (seq.step()) {
      seq.fillWithIndex(idx); 
      tmFluid.setPtr(ptr, idx);
      // only subtract from the diagonal parts. There should be a "correct" way to do this in 3d. FIXME
      src_in_out[idxr.getIndex(idx) + 0*vol] -= avgvT;
      src_in_out[idxr.getIndex(idx) + 3*vol] -= avgvT;
      src_in_out[idxr.getIndex(idx) + 5*vol] -= avgvT;
    }
    seq.reset();
    // compute forward transform
    fftw_execute(f_plan);

    double vt = std::sqrt(avgvT);
    double edt;
    // do relaxation    
    while (seq.step()) {
      seq.fillWithIndex(idx); 
      edt = std::exp(-kabs[idxr.getIndex(idx)]*vt*dt/2.0/Lucee::PI);
      for (unsigned i = 0; i < 6; ++i){
        /// hammett perkins says dP/dt = -n0 2 sqrt(2/pi) * k*vt T_k. Assume n0 is constant lah. 
        // so the Linear solution is T_new = T_old*exp(-k v factor dt)
        // then we need to add this back to the real solution
        // we are working in units of T/m, but since m is constant the solution is the same
        sol_in_out[idxr.getIndex(idx) + i*vol] = src_in_out[idxr.getIndex(idx) + i*vol]*edt;
      }
    }
    // compute inverse transform
    fftw_execute(b_plan);

    seq.reset();

    while (seq.step()) {
      seq.fillWithIndex(idx); 
      tmFluid.setPtr(ptr, idx);
      // only subtract from the diagonal parts. There should be a "correct" way to do this in 3d. FIXME
      double r = ptr[0];
      double u = ptr[1]/r;
      double v = ptr[2]/r;
      double w = ptr[3]/r;

      // update the pressure. P_ij = p_0 d_ij + n_0*T_ij_fluctuating + n m v v
      // division by volume due to FFTW not normalising output
      ptr[4] = avgvT*r + r*sol_in_out[idxr.getIndex(idx) + 0].real()/vol + r*u*u;
      ptr[5] = r*sol_in_out[idxr.getIndex(idx) + 1*vol].real()/vol + r*u*v;
      ptr[6] = r*sol_in_out[idxr.getIndex(idx) + 2*vol].real()/vol + r*u*w;
      ptr[7] = avgvT*r + r*sol_in_out[idxr.getIndex(idx) + 3*vol].real()/vol + r*v*v;
      ptr[8] = r*sol_in_out[idxr.getIndex(idx) + 4*vol].real()/vol + r*v*w;
      ptr[9] = avgvT*r + r*sol_in_out[idxr.getIndex(idx) + 5*vol].real()/vol + r*w*w;
    }
    seq.reset();

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>  
  void
  PeriodicCollisionlessHeatFluxUpdater<NDIM>::declareTypes()
  {
    //    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template class PeriodicCollisionlessHeatFluxUpdater<1>;
  template class PeriodicCollisionlessHeatFluxUpdater<2>;
  template class PeriodicCollisionlessHeatFluxUpdater<3>;

}

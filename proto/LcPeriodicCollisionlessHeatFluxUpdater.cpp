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
  static const unsigned PSIZE = 6; 

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

    int mpi_size = grid.getComm()->getNumProcs();
    int mpi_rank = grid.getComm()->getRank();
    std::cout<< "rank is " << mpi_rank << "\n"; // serial returns zero
    // for now disallow odd num procs
    if (mpi_rank%2 != 0 && mpi_rank != 1) {
      Lucee::Except lce("PeriodicCollisionlessHeatFluxUpdater needs even number of processors");
      throw lce;
    }

#ifdef HAVE_MPI
    ptrdiff_t N[NDIM]; // Global grid size for parallel version
    ptrdiff_t alloc_local, local_n0; 
    // here we only implement the 2-d version.
    for (unsigned i = 0; i < NDIM; ++i){
      N[i] = grid.getNumCells(i);
    }
    alloc_local = fftw_mpi_local_size_many(NDIM, N, 6, FFTW_MPI_DEFAULT_BLOCK,MPI_COMM_WORLD, &local_n0, &local_0_start);
#else
    int N[NDIM]; // Global grid size for parallel version
    for (unsigned i = 0; i < NDIM; ++i){
      N[i] = grid.getNumCells(i);
    }
    local_0_start = 0;
#endif    
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    int vol = grid.getGlobalRegion().getVolume();
    std::cout<< "volume is " << vol << "\n";
    
    /*    if (alloc_local != vol*6) {
      Lucee::Except lce("PeriodicCollisionlessHeatFluxUpdater has different memory requirements?!");
      throw lce;
      }*/
// construct arrays for use in FFTW -- no local allocation in the serial version
// using the FFTW many interface so allocate memory for all 6 arrays.
    src_in_out.resize(vol*6); // FIXMEP for parallel use alloc local.
    sol_in_out.resize(vol*6);


// computational space
    Lucee::Region<NDIM, double> compRgn = grid.getComputationalSpace();
    double Lx = compRgn.getShape(0);
    double Ly = compRgn.getShape(1);

// construct wave-number arrays. Funkyness is needed as DC component
// is stored in location [0]. FIXMEP for 3 dimensions.
    int nx = localRgn.getShape(0);
    kx.resize(nx);
    if (mpi_rank <= 1) {
      for (unsigned i=0; i<nx/2; ++i)
        kx[i] = 2*Lucee::PI/Lx*i;
      for (unsigned i=nx/2; i>0; --i)
        kx[nx-i] = -2*Lucee::PI/Lx*i;
    } else {
      if (2*mpi_rank < mpi_size) {
        for (unsigned i=0; i<nx; ++i)
          kx[i] = 2*Lucee::PI/Lx*(i+local_0_start);
      } else {
        for (unsigned i=0; i<nx; ++i)
        kx[i] = -2*Lucee::PI/Lx*(N[0] - local_0_start - i);
      }
    }
    // there is no decomposition in the y direction
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
    int howmany = PSIZE; // transform all 6 components
    int istride = howmany; // distance between elements in memory. For the parallel version we need stride = howmany.
    int ostride = howmany; 
    int idist = 1; // distance between first element of array -- interleaved data
  
    int odist = 1; 
#ifdef HAVE_MPI
  //  f_plan = fftw_mpi_plan_many_dft(NDIM, N, howmany, FFT_MPI_DEFAULT_BLOCK, FFT_MPI_DEFAULT_BLOCK, f_src_in_out, f_src_in_out, grid.getComm()->getMpiComm(), FFTW_FORWARD, FFTW_MEASURE);
  //  b_plan = fftw_mpi_plan_many_dft(NDIM, N, howmany, FFT_MPI_DEFAULT_BLOCK, FFT_MPI_DEFAULT_BLOCK, f_sol_in_out, f_sol_in_out, grid.getComm()->getMpiComm(), FFTW_FORWARD, FFTW_MEASURE);
#else
    f_plan = fftw_plan_many_dft(
             NDIM, N, howmany, f_src_in_out, NULL, istride, idist, f_src_in_out, NULL, ostride, odist, FFTW_FORWARD, FFTW_MEASURE); 
    b_plan = fftw_plan_many_dft(
             NDIM, N, howmany, f_sol_in_out, NULL, istride, idist, f_sol_in_out, NULL, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);
#endif

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
    int vol = grid.getGlobalRegion().getVolume();
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
      for (unsigned i = 0; i < PSIZE; ++i){
        // interleaved data
        src_in_out[idxr.getIndex(idx)*PSIZE + i] = vTarr[i];
      }
    }

    double avgvT, avgvTxx, avgvTyy, avgvTzz;
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &intIsovT, &avgvT, TX_SUM); 
    comm->allreduce(1, &intvTxx, &avgvTxx, TX_SUM); 
    comm->allreduce(1, &intvTyy, &avgvTyy, TX_SUM); 
    comm->allreduce(1, &intvTzz, &avgvTzz, TX_SUM); 
    avgvT = avgvT/vol; 
    avgvTxx = avgvTxx/vol; 
    avgvTyy = avgvTyy/vol; 
    avgvTzz = avgvTzz/vol; 
    // subtract isotopric temperature from input data
    seq.reset(); // reset the sequencer
      
    while (seq.step()) {
      seq.fillWithIndex(idx); 
      tmFluid.setPtr(ptr, idx);
      // only subtract from the diagonal parts. There should be a "correct" way to do this in 3d. FIXME
      // interleaved data
      src_in_out[idxr.getIndex(idx)*PSIZE ] -= avgvT;
      src_in_out[idxr.getIndex(idx)*PSIZE + 3] -= avgvT;
      src_in_out[idxr.getIndex(idx)*PSIZE + 5] -= avgvT;
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
      // should there be a factor of 2*PI here?
      for (unsigned i = 0; i < PSIZE; ++i){
        /// hammett perkins says dP/dt = -n0 2 sqrt(2/pi) * k*vt T_k. Assume n0 is constant lah. 
        // so the Linear solution is T_new = T_old*exp(-k v factor dt)
        // then we need to add this back to the real solution
        // we are working in units of T/m, but since m is constant the solution is the same
        sol_in_out[idxr.getIndex(idx)*PSIZE + i] = src_in_out[idxr.getIndex(idx)*PSIZE + i]*edt;
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
      ptr[4] = avgvT*r + r*sol_in_out[idxr.getIndex(idx)*PSIZE].real()/vol + r*u*u;
      ptr[5] = r*sol_in_out[idxr.getIndex(idx)*PSIZE + 1].real()/vol + r*u*v;
      ptr[6] = r*sol_in_out[idxr.getIndex(idx)*PSIZE + 2].real()/vol + r*u*w;
      ptr[7] = avgvT*r + r*sol_in_out[idxr.getIndex(idx)*PSIZE + 3].real()/vol + r*v*v;
      ptr[8] = r*sol_in_out[idxr.getIndex(idx)*PSIZE + 4].real()/vol + r*v*w;
      ptr[9] = avgvT*r + r*sol_in_out[idxr.getIndex(idx)*PSIZE + 5].real()/vol + r*w*w;
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

/**
 * @file	LcPeriodicPoisson2DUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with periodic BCs
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMathPhysConstants.h>
#include <LcPeriodicParallelPoisson2DUpdater.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

namespace Lucee
{

/** Class id: this is used by registration system */
  const char *PeriodicParallelPoisson2DUpdater::id = "PeriodicParallelPoisson2D";

  PeriodicParallelPoisson2DUpdater::PeriodicParallelPoisson2DUpdater()
    : UpdaterIfc()
  {
  }

  void
  PeriodicParallelPoisson2DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
  }

  void
  PeriodicParallelPoisson2DUpdater::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
    // This is a first attempt to use FFTW's parallel functionality to calculate the transform of the 2d poisson equation. Because FFTW does the decomposition by itself it is necessary to give it the size of the entire grid. The assumption is that data are decomposed in slabs
    // For now we assume that the two decompositions are the same. 
// get hold of grid
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();
    // global size needed by FFTW
    const ptrdiff_t NX = grid.getNumCells(0);
    const ptrdiff_t NY = grid.getNumCells(1); 
    // fftw variables for local use
    ptrdiff_t alloc_local, local_n0;
    // move this to gkyell.cxx later.
    fftw_mpi_init();
    // create the fft local data size (hopefully grid.getComm() actually does something.
    alloc_local = fftw_mpi_local_size_2d(NX , NY, MPI_COMM_WORLD, &local_n0, &local_0_start);
    // FFTW may want more tha nx*ny 
    src_in_out.resize(alloc_local);
    sol_in_out.resize(alloc_local);

// local region to index
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
    int vol = localRgn.getVolume();

    
// computational space
    Lucee::Region<2, double> compRgn = grid.getComputationalSpace();
    double Lx = compRgn.getShape(0);
    double Ly = compRgn.getShape(1);

// construct wave-number arrays. Funkyness is needed as DC component
// is stored in location [0].
// In the parallel version we also need to get the correct wavenumbers in each process. Cuts are in the x direction
    int nx = localRgn.getShape(0);
    int mpi_size = grid.getComm()->getNumProcs();
    int mpi_rank = grid.getComm()->getRank();

    kx.resize(nx);
    if (2*mpi_rank < mpi_size) {
      for (unsigned i=0; i<nx; ++i)
        kx[i] = 2*Lucee::PI/Lx*(i+local_0_start);
      if (local_0_start == 0){
        kx[0] = 1e-8;
      }
    } else {
      for (unsigned i=0; i<nx; ++i)
        kx[i] = -2*Lucee::PI/Lx*(NX - local_0_start - i);
    }

    // y decomposition is the same
    int ny = localRgn.getShape(1);
    ky.resize(ny);
    for (unsigned i=0; i<ny/2; ++i)
      ky[i] = 2*Lucee::PI/Ly*i;
    for (unsigned i=ny/2; i>0; --i)
      ky[ny-i] = -2*Lucee::PI/Ly*i;
    // reset DC component to avoid divide-by-zero
    ky[0] = 1.0e-8;  

// these ugly casts are needed so we can pass pointers of correct type
// to FFTW. This works because for *most* C++ compilers
// complex<double> is the same as double [2] which FFTW uses and is
// typedef-ed as fftw_complex. It is possible that this can cause
// problems if some C++ compiler does not support this.
    f_src_in_out = reinterpret_cast<fftw_complex*> (&src_in_out[0]);
    f_sol_in_out = reinterpret_cast<fftw_complex*> (&sol_in_out[0]);

// now create plans for forward and inverse FFTs
    f_plan = fftw_mpi_plan_dft_2d(NX, NY, f_src_in_out, f_src_in_out, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
    b_plan = fftw_mpi_plan_dft_2d(NX, NY, f_sol_in_out, f_sol_in_out, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  Lucee::UpdaterStatus
  PeriodicParallelPoisson2DUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input/output fields
    const Lucee::Field<2, double>& src = this->getInp<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& sol = this->getOut<Lucee::Field<2, double> >(0);

// local indexing region
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
    int vol = localRgn.getVolume();

// compute integrated source to remove from source. This is needed to
// ensure a properly posed problem as we are working on a periodic
// domain.
    double intSrc = 0.0;
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
        intSrc += src(i,j,0);

    double globalSum;
    int globalVol;
    MPI_Allreduce(&intSrc, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vol, &globalVol, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    std::cout<< globalSum/globalVol << "\n";
// copy source into FFTW in/out array
    Lucee::RowMajorIndexer<2> idxr(localRgn);
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i){
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j){
        src_in_out[idxr.getIndex(i,j)] = src(i,j,0) - globalSum/globalVol;
        //        f_src_in_out[idxr.getIndex(i,j)][1] = 0;
      }
    }

// compute forward transform
    fftw_execute(f_plan);
// now compute FFT of solution
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j){
        sol_in_out[idxr.getIndex(i,j)] 
          = -src_in_out[idxr.getIndex(i,j)]/(kx[i-local_0_start]*kx[i-local_0_start] + ky[j]*ky[j]);
        //        f_sol_in_out[idxr.getIndex(i,j)][1] 
        //          = -f_src_in_out[idxr.getIndex(i,j)][1]/(kx[i-local_0_start]*kx[i-local_0_start] + ky[j]*ky[j]);
      }
// compute inverse transform
    fftw_execute(b_plan);

// store solution (division by volume is needed as the FFTs are not normalized)
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
        sol(i,j,0) = sol_in_out[idxr.getIndex(i,j)].real()/globalVol;

// apply periodic BCs to solution array
    sol.applyPeriodicBc(0); 
    sol.applyPeriodicBc(1);

    return Lucee::UpdaterStatus();
  }
  
  void
  PeriodicParallelPoisson2DUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}

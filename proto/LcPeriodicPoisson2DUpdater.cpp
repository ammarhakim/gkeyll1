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
#include <LcPeriodicPoisson2DUpdater.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

namespace Lucee
{

/** Class id: this is used by registration system */
  const char *PeriodicPoisson2DUpdater::id = "PeriodicPoisson2D";

  PeriodicPoisson2DUpdater::PeriodicPoisson2DUpdater()
    : UpdaterIfc()
  {
  }

  void
  PeriodicPoisson2DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
  }

  void
  PeriodicPoisson2DUpdater::initialize()
  {
// call base class method
    UpdaterIfc::initialize();

// get hold of grid
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();
// local region to index
    Lucee::Region<2, int> localRgn = grid.getLocalBox();
    int vol = localRgn.getVolume();

// computational space
    Lucee::Region<2, double> compRgn = grid.getComputationalSpace();
    double Lx = compRgn.getShape(0);
    double Ly = compRgn.getShape(1);

// construct arrays for use in FFTW
    src_in_out.resize(vol);
    sol_in_out.resize(vol);

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

// reset DC component to avoid divide-by-zero
    kx[0] = 1.0e-8;
    ky[0] = 1.0e-8;

// these ugly casts are needed so we can pass pointers of correct type
// to FFTW. This works because for *most* C++ compilers
// complex<double> is the same as double [2] which FFTW uses and is
// typedef-ed as fftw_complex. It is possible that this can cause
// problems if some C++ compiler does not support this.
    f_src_in_out = reinterpret_cast<fftw_complex*> (&src_in_out[0]);
    f_sol_in_out = reinterpret_cast<fftw_complex*> (&sol_in_out[0]);

// now create plans for forward and inverse FFTs
    f_plan = fftw_plan_dft_2d(
      nx, ny, f_src_in_out, f_src_in_out, FFTW_FORWARD, FFTW_ESTIMATE);
    b_plan = fftw_plan_dft_2d(
      nx, ny, f_sol_in_out, f_sol_in_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  Lucee::UpdaterStatus
  PeriodicPoisson2DUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input/output fields
    const Lucee::Field<2, double>& src = this->getInp<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& sol = this->getOut<Lucee::Field<2, double> >(0);

// local indexing region
    Lucee::Region<2, int> localRgn = grid.getLocalBox();
    int vol = localRgn.getVolume();

// compute integrated source to remove from source. This is needed to
// ensure a properly posed problem as we are working on a periodic
// domain.
    double intSrc = 0.0;
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
        intSrc += src(i,j,0);

// copy source into FFTW in/out array
    Lucee::RowMajorIndexer<2> idxr(localRgn);
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
        src_in_out[idxr.getIndex(i,j)] = src(i,j,0) - intSrc/vol;

// compute forward transform
    fftw_execute(f_plan);

// now compute FFT of solution
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
        sol_in_out[idxr.getIndex(i,j)] 
          = -src_in_out[idxr.getIndex(i,j)]/(kx[i]*kx[i] + ky[j]*ky[j]);

// compute inverse transform
    fftw_execute(b_plan);

// store solution (division by volume is needed as the FFTs are not normalized)
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
        sol(i,j,0) = sol_in_out[idxr.getIndex(i,j)].real()/vol;

// apply periodic BCs to solution array
    sol.applyPeriodicBc(0); 
    sol.applyPeriodicBc(1);

    return Lucee::UpdaterStatus();
  }
  
  void
  PeriodicPoisson2DUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}

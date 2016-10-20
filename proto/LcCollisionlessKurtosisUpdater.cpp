/**
 * @file	LcCollisionlessKurtosisUpdater.cpp
 *
 * @brief	Updater for relaxation of the nonlocal heat flux tensor
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// txbase includes
#include <TxCommBase.h>
#ifdef HAVE_MPI
#include <TxMpiBase.h>
#endif

// lucee includes
#include <LcMathPhysConstants.h>
#include <LcCollisionlessKurtosisUpdater.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

// Eigen includes
#include <Eigen/Eigen>


namespace Lucee
{
  static const unsigned PSIZE = 10; //ten tensor components
  static const unsigned OFF = 10; // first heat flux component is at location 10
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

/** Class id: this is used by registration system */
  template <> const char *CollisionlessKurtosisUpdater<1>::id = "CollisionlessKurtosis1D";
  template <> const char *CollisionlessKurtosisUpdater<2>::id = "CollisionlessKurtosis2D";
  template <> const char *CollisionlessKurtosisUpdater<3>::id = "CollisionlessKurtosis3D";

  template <unsigned NDIM>
  CollisionlessKurtosisUpdater<NDIM>::CollisionlessKurtosisUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  CollisionlessKurtosisUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
    if (tbl.hasNumber("factor")){
      scale_factor = tbl.getNumber("factor");
    } else {
      scale_factor = 1.0;
    }
    if (tbl.hasString("closure")){
      std::string cl = tbl.getString("closure");
      if (cl == "cai"){
        closure = CL_CAI;
      } else if (cl == "conservative"){
        closure = CL_CONSERVATIVE;
      } else {
        Lucee::Except lce("CollisionlessKurtosisUpdater::readInput: 'closure' ");
        lce << cl << " not recognized!" << std::endl;
        throw lce;
      }
    } else {
      closure = CL_CAI;
    }

    // determine the solver method
    if (tbl.hasString("solver")){
      std::string sl = tbl.getString("solver");
      if (sl == "ndim"){
        solver_method = S_ND;
      } else {
        solver_method = S_NAIVE;
      }
    } else {
      solver_method = S_NAIVE;
    }
  }

  template <unsigned NDIM>
  void
  CollisionlessKurtosisUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
// get hold of grid
    const Lucee::RectCartGrid<NDIM>& grid 
      = this->getGrid<Lucee::RectCartGrid<NDIM> >();

// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    int vol = grid.getGlobalRegion().getVolume();
    std::cout<< "volume is " << vol << "\n";

    // temporary restriction of number of processors due to interfacing with FFTW.
    int mpi_size = grid.getComm()->getNumProcs();
    int mpi_rank = grid.getComm()->getRank();
    std::cout<< "rank is " << mpi_rank << "\n"; // serial returns zero
    // for now disallow odd num procs
    if (mpi_size%2 != 0 && mpi_size != 1) {
      Lucee::Except lce("CollisionlessKurtosisUpdater needs even number of processors");
      throw lce;
    }

#ifdef HAVE_MPI
    TxMpiBase *comm = static_cast<TxMpiBase*>(this->getComm());
    ptrdiff_t N[NDIM]; // Global grid size for parallel version
    ptrdiff_t alloc_local, local_n0; 
    // here we only implement the 2-d version.
    for (unsigned i = 0; i < NDIM; ++i){
      N[i] = grid.getNumCells(i);
    }
    alloc_local = fftw_mpi_local_size_many(NDIM, N, PSIZE, FFTW_MPI_DEFAULT_BLOCK,comm->getMpiComm(), &local_n0, &local_0_start);
    // the parallel allocation can be larger than the volume*10
    src_in_out.resize(alloc_local); 
    sol_in_out.resize(alloc_local);
    

#else
    int N[NDIM]; // Global grid size for parallel version
    for (unsigned i = 0; i < NDIM; ++i){
      N[i] = grid.getNumCells(i);
    }
    local_0_start = 0;
    src_in_out.resize(vol*10);
    sol_in_out.resize(vol*10);

#endif    
    
// construct arrays for use in FFTW -- no local allocation in the serial version
// using the FFTW many interface so allocate memory for all 10 arrays.


// computational space
    Lucee::Region<NDIM, double> compRgn = grid.getComputationalSpace();
    double Lx = compRgn.getShape(0);

    int nLoc[3];
    double lGlobal[3]; // global physical dimensions
    for (unsigned i = 0; i < NDIM; ++i){
      nLoc[i] = localRgn.getShape(i);
      lGlobal[i] = compRgn.getShape(i);
    } 
    for (unsigned i = NDIM; i < 3; ++i){
      nLoc[i] = 1; // fill in the rest of the dims with 1
      lGlobal[i] = 1.0; //doesn't matter but not zero
    }

    int nx = nLoc[0];
    kx.resize(nx);
    for (unsigned i=0; i<nx; ++i){
      if ((i+local_0_start) <= N[0]/2){
        // positive wavenumber
        kx[i] = 2*Lucee::PI/Lx*(i+local_0_start);
      } else { 
        // negative wavenumber
        kx[i] = -2*Lucee::PI/Lx*(N[0] - local_0_start - i);
      }
    }
    //    if (local_0_start == 0){
      //      kx[0] = 1e-8; // avoid division by zero
    //    }
    
    // there is no decomposition in the y direction
    double Ly = lGlobal[1];
    int ny = nLoc[1];
    ky.resize(ny);
    for (unsigned i=0; i<ny; ++i){
      if (i<= ny/2){
        ky[i] = 2*Lucee::PI/Ly*i;
      } else {
        ky[i] = -2*Lucee::PI/Ly*(ny - i);
      }
    }

    double Lz = lGlobal[2];
    int nz = nLoc[2];
    kz.resize(nz);
    for (unsigned i=0; i < nz; ++i) {
      if (i <= nz/2) {
        kz[i] = 2*Lucee::PI/Lz*i;
      } else {
        kz[i] = -2*Lucee::PI/Lz*(nz-i);
      }
    }
    //    for (int i=nz/2-1; i > 0; --i)
    //      kz[nz-i] = -2*Lucee::PI/Lz*i;
    //    kz[0] = 1e-8; // avoid division by zero

    // store |k| for collisionless relaxation.
    kabs.resize(nx*ny*nz); 
    for (unsigned i=0; i<nx; ++i){
      for (unsigned j=0; j<ny; ++j){
        for (unsigned k=0; k<nz; ++k) {
          kabs[k + nz*(j+ny*i)] = std::sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
        }
      }
    }
    // prevent division by zero!
    if (local_0_start == 0){
      kabs[0] = 2*Lucee::PI/Lx*1e-3;
    }


// these ugly casts are needed so we can pass pointers of correct type
// to FFTW. This works because for *most* C++ compilers
// complex<double> is the same as double [2] which FFTW uses and is
// typedef-ed as fftw_complex. It is possible that this can cause
// problems if some C++ compiler does not support this.
    f_src_in_out = reinterpret_cast<fftw_complex*> (&src_in_out[0]);
    f_sol_in_out = reinterpret_cast<fftw_complex*> (&sol_in_out[0]);
  
// now create plans for forward and inverse FFTs
// fft_plan_many_dft(int rank, const int *n, int howmany, fftw_complex *in, const int *inembed, int istride, int idist, fftw_complex * out, const int *oembed, int ostride, int odist, int sign, unsigned flags); 
    int howmany = PSIZE; // transform all 10 components
    int istride = howmany; // distance between elements in memory. For the parallel version we need stride = howmany. Thus the serial version has the same format.
    int ostride = howmany; 
    int idist = 1; // distance between first element of each array in memory-- interleaved data
    int odist = 1; 
#ifdef HAVE_MPI
    // FFTW_MEASURE gives better speed at the cost of setup time. FFTW_PATIENT is supposed to be even better than MEASURE, but it does not seem to have an effect for small systems.
    f_plan = fftw_mpi_plan_many_dft(NDIM, N, howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, f_src_in_out, f_src_in_out, comm->getMpiComm(), FFTW_FORWARD, FFTW_MEASURE);
    b_plan = fftw_mpi_plan_many_dft(NDIM, N, howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, f_sol_in_out, f_sol_in_out, comm->getMpiComm(), FFTW_BACKWARD, FFTW_MEASURE);
#else
    f_plan = fftw_plan_many_dft(
             NDIM, N, howmany, f_src_in_out, NULL, istride, idist, f_src_in_out, NULL, ostride, odist, FFTW_FORWARD, FFTW_MEASURE); 
    b_plan = fftw_plan_many_dft(
             NDIM, N, howmany, f_sol_in_out, NULL, istride, idist, f_sol_in_out, NULL, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);
#endif

  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  CollisionlessKurtosisUpdater<NDIM>::update(double t)
  {
// get hold of grid and find dt
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double dt = t-this->getCurrTime();
// get fields
    Lucee::Field<NDIM, double>& tmFluid = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> ptr = tmFluid.createPtr();
    int idx[NDIM]; 
// local indexing region
    Lucee::Region<NDIM, int> localRgn = tmFluid.getRegion();
    Lucee::RowMajorSequencer <NDIM> seq(localRgn);
    Lucee::RowMajorIndexer<NDIM> idxr(localRgn);
    int vol = grid.getGlobalRegion().getVolume();
    double qarr[PSIZE];

    double intIsovT = 0.0;
    // do stepping
    while (seq.step()) { 

      seq.fillWithIndex(idx);
      tmFluid.setPtr(ptr, idx); 
      // we need this for calculating thermal velocity
      double r = ptr[0];
      double u = ptr[1]/r;
      double v = ptr[2]/r;
      double w = ptr[3]/r;
      double vTxx = (ptr[4]-r*u*u)/r;
      double vTyy = (ptr[7]-r*v*v)/r;
      double vTzz = (ptr[9]-r*w*w)/r;
      if (closure == CL_CONSERVATIVE){
        double vTxy = (ptr[5]-r*u*v)/r;
        double vTxz = (ptr[6]-r*u*w)/r;
        double vTyz = (ptr[8]-r*v*w)/r;
        qarr[0] = ptr[Q111] - 3*r*vTxx*u - r*u*u*u;
        qarr[1] = ptr[Q112] - 2*r*vTxy*u - r*vTxx*v - r*u*u*v;
        qarr[2] = ptr[Q113] - 2*r*vTxz*u - r*vTxx*w - r*u*u*w;
        qarr[3] = ptr[Q122] - r*vTyy*u - 2*r*vTxy*v - r*u*v*v;
        qarr[4] = ptr[Q123] - r*vTyz*u - r*vTxz*v - r*vTxy*w - r*u*v*w;
        qarr[5] = ptr[Q133] - r*vTzz*u - 2*r*vTxz*w - r*u*w*w;
        qarr[6] = ptr[Q222] - 3*r*vTyy*v - r*v*v*v;
        qarr[7] = ptr[Q223] - 2*r*vTyz*v - r*vTyy*w - r*v*v*w;
        qarr[8] = ptr[Q233] - r*vTzz*v - 2*r*vTyz*w - r*v*w*w;
        qarr[9] = ptr[Q333] - 3*r*vTzz*w - r*w*w*w;
      }
      // find global isotropic temperature and fill matrix
      intIsovT += (vTxx+vTyy+vTzz)/3.0;
      // fill the matrix for transformation
      for (unsigned i = 0; i < PSIZE; ++i){
        // interleaved data

        if (closure == CL_CONSERVATIVE){
          src_in_out[idxr.getIndex(idx)*PSIZE + i] = qarr[i];
        } else {
          src_in_out[idxr.getIndex(idx)*PSIZE + i] = ptr[i+OFF];
        }
      }
    }
    seq.reset(); // reset the sequencer
    // calculate the thermal speed
    double avgvT;
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &intIsovT, &avgvT, TX_SUM); 
    avgvT = avgvT/vol; 
    // compute forward transform
    fftw_execute(f_plan);

    double vt = std::sqrt(avgvT);
    double edt;

    /* eigen matrices for point implicit time step */
    Eigen::MatrixXcd prLhs;
    prLhs = Eigen::MatrixXcd::Constant(10, 10, 0.0);
    Eigen::VectorXcd prRhs(10);
    Eigen::VectorXcd prSol(10);

    // do relaxation proportional to wavenumber
    while (seq.step()) {
      seq.fillWithIndex(idx); 
      // in the 1-D hammett-perkins dr has a component proportional to T as well. Here we ignore it and use dq/dt = |k| D q
      if (solver_method != S_ND){ 
        edt = std::exp(-kabs[idxr.getIndex(idx)]*vt*dt*scale_factor*std::sqrt(8*Lucee::PI)/(3*Lucee::PI - 8.0));
        for (unsigned i = 0; i < PSIZE; ++i){
          sol_in_out[idxr.getIndex(idx)*PSIZE + i] = src_in_out[idxr.getIndex(idx)*PSIZE + i]*edt;
        }
      } else { // the 3d matrix method
        for (unsigned i=0; i < PSIZE; ++i){
          prRhs(i) = src_in_out[idxr.getIndex(idx)*PSIZE+i];
        }

        double kxLoc = kx[idx[0]-local_0_start];
        double kyLoc = ky[idx[1]];
        double kzLoc;
        double kLoc = kabs[idxr.getIndex(idx)];
        double chivtdt = std::sqrt(8*Lucee::PI)/(3*Lucee::PI-8.0)/kLoc*vt*dt*kLoc*scale_factor; // in general the k vt should be replaced by sqrt(ki tij kj)
        if (NDIM < 3) { // there should be no reason to call this for N < 2
          kzLoc = 0.0;
        } else {
          kzLoc = kz[idx[2]];
        }
        // coeffiecient matrix for implicit solve
        prLhs(0,0) = 1 + (chivtdt*((kxLoc*kxLoc) + (kyLoc*kyLoc)/4 + (kzLoc*kzLoc)/4))/kLoc; prLhs(0,1) = (3*chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(0,2) = (3*chivtdt*kxLoc*kzLoc)/(4*kLoc); 
        prLhs(1,1) = 1 + (chivtdt*((3*(kxLoc*kxLoc))/4 + (kyLoc*kyLoc)/2 + (kzLoc*kzLoc)/4))/kLoc; prLhs(1,0) = (chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(1,3) = (chivtdt*kxLoc*kyLoc)/(2*kLoc); prLhs(1,2) = (chivtdt*kyLoc*kzLoc)/(4*kLoc); prLhs(1,4) = (chivtdt*kxLoc*kzLoc)/(2*kLoc); 
        prLhs(2,2) = 1 + (chivtdt*((3*(kxLoc*kxLoc))/4 + (kyLoc*kyLoc)/4 + (kzLoc*kzLoc)/2))/kLoc; prLhs(2,0) = (chivtdt*kxLoc*kzLoc)/(4*kLoc); prLhs(2,1) = (chivtdt*kyLoc*kzLoc)/(4*kLoc); prLhs(2,4) = (chivtdt*kxLoc*kyLoc)/(2*kLoc); prLhs(2,5) = (chivtdt*kxLoc*kzLoc)/(2*kLoc); 
        prLhs(3,3) = 1 + (chivtdt*((kxLoc*kxLoc)/2 + (3*(kyLoc*kyLoc))/4 + (kzLoc*kzLoc)/4))/kLoc; prLhs(3,1) = (chivtdt*kxLoc*kyLoc)/(2*kLoc); prLhs(3,6) = (chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(3,4) = (chivtdt*kyLoc*kzLoc)/(2*kLoc); prLhs(3,7) = (chivtdt*kxLoc*kzLoc)/(4*kLoc); 
        prLhs(4,4) = 1 + (chivtdt*((kxLoc*kxLoc)/2 + (kyLoc*kyLoc)/2 + (kzLoc*kzLoc)/2))/kLoc; prLhs(4,1) = (chivtdt*kxLoc*kzLoc)/(4*kLoc); prLhs(4,2) = (chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(4,3) = (chivtdt*kyLoc*kzLoc)/(4*kLoc); prLhs(4,7) = (chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(4,5) = (chivtdt*kyLoc*kzLoc)/(4*kLoc); prLhs(4,8) = (chivtdt*kxLoc*kzLoc)/(4*kLoc); 
        prLhs(5,5) = 1 + (chivtdt*((kxLoc*kxLoc)/2 + (kyLoc*kyLoc)/4 + (3*(kzLoc*kzLoc))/4))/kLoc; prLhs(5,2) = (chivtdt*kxLoc*kzLoc)/(2*kLoc); prLhs(5,4) = (chivtdt*kyLoc*kzLoc)/(2*kLoc); prLhs(5,8) = (chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(5,9) = (chivtdt*kxLoc*kzLoc)/(4*kLoc); 
        prLhs(6,6) = 1 + (chivtdt*((kxLoc*kxLoc)/4 + (kyLoc*kyLoc) + (kzLoc*kzLoc)/4))/kLoc; prLhs(6,3) = (3*chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(6,7) = (3*chivtdt*kyLoc*kzLoc)/(4*kLoc); 
        prLhs(7,7) = 1 + (chivtdt*((kxLoc*kxLoc)/4 + (3*(kyLoc*kyLoc))/4 + (kzLoc*kzLoc)/2))/kLoc; prLhs(7,3) = (chivtdt*kxLoc*kzLoc)/(4*kLoc); prLhs(7,4) = (chivtdt*kxLoc*kyLoc)/(2*kLoc); prLhs(7,6) = (chivtdt*kyLoc*kzLoc)/(4*kLoc); prLhs(7,8) = (chivtdt*kyLoc*kzLoc)/(2*kLoc); 
        prLhs(8,8) = 1 + (chivtdt*((kxLoc*kxLoc)/4 + (kyLoc*kyLoc)/2 + (3*(kzLoc*kzLoc))/4))/kLoc; prLhs(8,4) = (chivtdt*kxLoc*kzLoc)/(2*kLoc); prLhs(8,5) = (chivtdt*kxLoc*kyLoc)/(4*kLoc); prLhs(8,7) = (chivtdt*kyLoc*kzLoc)/(2*kLoc); prLhs(8,9) = (chivtdt*kyLoc*kzLoc)/(4*kLoc); 
        prLhs(9,9) = 1 + (chivtdt*((kxLoc*kxLoc)/4 + (kyLoc*kyLoc)/4 + (kzLoc*kzLoc)))/kLoc; prLhs(9,5) = (3*chivtdt*kxLoc*kzLoc)/(4*kLoc); prLhs(9,8) = (3*chivtdt*kyLoc*kzLoc)/(4*kLoc); 
        prSol = prLhs.partialPivLu().solve(prRhs);
        for (unsigned i = 0; i < PSIZE; ++i){
          sol_in_out[idxr.getIndex(idx)*PSIZE + i] = prSol(i);
        }
      }
    }
    // compute inverse transform
    fftw_execute(b_plan);

    seq.reset();

    // update the pressure by adding the damped fluctuations
    while (seq.step()) {
      seq.fillWithIndex(idx); 
      tmFluid.setPtr(ptr, idx);
      double qshift[PSIZE];
      if (closure == CL_CONSERVATIVE){
        double r = ptr[0];
        double u = ptr[1]/r;
        double v = ptr[2]/r;
        double w = ptr[3]/r;
        double vTxx = (ptr[4]-r*u*u)/r;
        double vTyy = (ptr[7]-r*v*v)/r;
        double vTzz = (ptr[9]-r*w*w)/r;
        double vTxy = (ptr[5]-r*u*v)/r;
        double vTxz = (ptr[6]-r*u*w)/r;
        double vTyz = (ptr[8]-r*v*w)/r;
        qshift[0] = -3*r*vTxx*u - r*u*u*u;
        qshift[1] = -2*r*vTxy*u - r*vTxx*v - r*u*u*v;
        qshift[2] = -2*r*vTxz*u - r*vTxx*w - r*u*u*w;
        qshift[3] = -r*vTyy*u - 2*r*vTxy*v - r*u*v*v;
        qshift[4] = -r*vTyz*u - r*vTxz*v - r*vTxy*w - r*u*v*w;
        qshift[5] = -r*vTzz*u - 2*r*vTxz*w - r*u*w*w;
        qshift[6] = -3*r*vTyy*v - r*v*v*v;
        qshift[7] = -2*r*vTyz*v - r*vTyy*w - r*v*v*w;
        qshift[8] = -r*vTzz*v - 2*r*vTyz*w - r*v*w*w;
        qshift[9] = -3*r*vTzz*w - r*w*w*w;
      }

      // update the pressure. P_ij = p_0 d_ij + n_0*T_ij_fluctuating + n m v v
      // division by volume due to FFTW not normalising output
      for (unsigned i = 0; i < PSIZE; i++){
        if (closure == CL_CONSERVATIVE) {
          ptr[i+OFF] = sol_in_out[idxr.getIndex(idx)*PSIZE + i].real()/vol - qshift[i];
        } else {
          ptr[i+OFF] = sol_in_out[idxr.getIndex(idx)*PSIZE + i].real()/vol;
        }
      }
    }
    seq.reset();

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>  
  void
  CollisionlessKurtosisUpdater<NDIM>::declareTypes()
  {
    //    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template class CollisionlessKurtosisUpdater<1>;
  template class CollisionlessKurtosisUpdater<2>;
  template class CollisionlessKurtosisUpdater<3>;

}

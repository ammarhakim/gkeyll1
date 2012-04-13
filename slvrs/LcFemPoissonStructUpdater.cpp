/**
 * @file	LcFemPoissonStructUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with FEM scheme on a structured grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFemPoissonStructUpdater.h>
#include <LcGlobals.h>
#include <LcMatrix.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

// txbase includes
#include <TxCommBase.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <sstream>
#include <vector>

namespace Lucee
{
  static const unsigned DX = 0;
  static const unsigned DY = 1;
  static const unsigned DZ = 2;
  static const unsigned LO = 0;
  static const unsigned HI = 1;
  static const unsigned NO_BC = 0;
  static const unsigned DIRICHLET_BC = 1;
  static const unsigned NEUMANN_BC = 2;

  template <> const char *FemPoissonStructUpdater<1>::id = "FemPoisson1D";
  template <> const char *FemPoissonStructUpdater<2>::id = "FemPoisson2D";
  template <> const char *FemPoissonStructUpdater<3>::id = "FemPoisson3D";

  template <unsigned NDIM>
  FemPoissonStructUpdater<NDIM>::FemPoissonStructUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  FemPoissonStructUpdater<NDIM>::~FemPoissonStructUpdater()
  { // get rid of stuff
    MatDestroy(stiffMat);
    VecDestroy(globalSrc);
    VecDestroy(initGuess);
  }

  template <unsigned NDIM>
  void
  FemPoissonStructUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except(
        "FemPoissonStructUpdater::readInput: Must specify element to use using 'basis'");

// check if source nodes are shared
    srcNodesShared = true;
    if (tbl.hasBool("sourceNodesShared"))
      srcNodesShared = tbl.getBool("sourceNodesShared");

    for (unsigned i=0; i<NDIM; ++i) 
      periodicFlgs[i] = false;

// check if any direction is periodic
    if (tbl.hasNumVec("periodicDirs"))
    {
      std::vector<double> pd = tbl.getNumVec("periodicDirs");
      for (unsigned i=0; i<pd.size(); ++i)
      {
        if ((pd[i]<0) || (pd[i]>=NDIM))
        {
          throw Lucee::Except("FemPoissonStructUpdater:readInput: Incorrect 'periodicDirs'");
        }
        periodicFlgs[(unsigned) pd[i]] = true;
      }
    }

// check if all directions are periodic
    allPeriodic = true;
    for (unsigned i=0; i<NDIM; ++i)
      allPeriodic = allPeriodic && periodicFlgs[i];

// get BCs to apply
    if (tbl.hasTable("bcLeft"))
    {
      bc[DX][LO] = getBcData(tbl.getTable("bcLeft"));
    }
    if (tbl.hasTable("bcRight"))
    {
      bc[DX][HI] = getBcData(tbl.getTable("bcRight"));
    }
    if (NDIM>1)
    {
      if (tbl.hasTable("bcBottom"))
      {
        bc[DY][LO] = getBcData(tbl.getTable("bcBottom"));
      }
      if (tbl.hasTable("bcTop"))
      {
        bc[DY][HI] = getBcData(tbl.getTable("bcTop"));
      }
    }
    if (NDIM>2)
    {
      if (tbl.hasTable("bcBack"))
      {
        bc[DZ][LO] = getBcData(tbl.getTable("bcBack"));
      }
      if (tbl.hasTable("bcFront"))
      {
        bc[DZ][HI] = getBcData(tbl.getTable("bcFront"));
      }
    }

// some sanity checks: must either specify both BCs along a direction
// or none. In the latter case the BC is assumed to be periodic.

    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      if (periodicFlgs[dir] == false)
      {
// ensure both BCs are specified
        if (!bc[dir][0].isSet || !bc[dir][1].isSet)
        {
          Lucee::Except lce("FemPoissonStructUpdater:: Must specify BCs on each side");
          throw lce;
        }
      }
      else
      {
// ensure no BCs are specified
        if (bc[dir][0].isSet || bc[dir][1].isSet)
        {
          Lucee::Except lce("FemPoissonStructUpdater:: Cannot specify BCs");
          lce << " if a direction is periodic";
          throw lce;
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  FemPoissonStructUpdater<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// number of global nodes
    unsigned nglobal = nodalBasis->getNumGlobalNodes();
// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// map to hold periodicially identified nodes
    std::map<int, int> periodicNodeMap;

#ifdef HAVE_MPI
    throw Lucee::Except("FemPoissonStructUpdater does not yet work in parallel!");
#else
// Explicit initialization of stiffness matrix speeds up
// initialization tremendously.
    int nz = 10; // number of non-zero entries per row (WHAT SHOULD IT REALLY BE?)
    if (nodalBasis->getNumNodes() > 4)
      nz = 21; // HACK FOR NOW
    MatCreateSeqAIJ(PETSC_COMM_SELF, nglobal, nglobal, nz, PETSC_NULL, &stiffMat);
#endif
    MatSetFromOptions(stiffMat);

// create and setup vector (NOTE: ACCORDING TO MIKE MCCOURT ONE SHOULD
// USE MatGetVecs INSTEAD OF VecCreate. THIS ENSURES THAT THE PARALLEL
// LAYOUT OF THE MATRIX & VECTOR ARE THE SAME)
    VecCreate(MPI_COMM_WORLD, &globalSrc);
    VecSetSizes(globalSrc, nglobal, PETSC_DECIDE);
    VecSetFromOptions(globalSrc);

// local stiffness matrix
    Lucee::Matrix<double> localStiff(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// storage for passing to petsc
    std::vector<PetscScalar> vals(nlocal*nlocal);

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// global region for grid
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
// create sequencer for looping over local box
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion(); 
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    int idx[NDIM];

// loop, creating stiffness matrix
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// set index into element basis
      nodalBasis->setIndex(idx);

// get local stiffness matrix
      nodalBasis->getStiffnessMatrix(localStiff);
// construct arrays for passing into Petsc
      for (unsigned k=0; k<nlocal; ++k)
      {
        for (unsigned m=0; m<nlocal; ++m)
          vals[nlocal*k+m] = -localStiff(k,m); // Default PetSc layout is row-major
      }

// get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

// insert into global stiffness matrix, adding them to existing value
      MatSetValues(stiffMat, nlocal, &lgMap[0], nlocal, &lgMap[0],
        &vals[0], ADD_VALUES);
    }

    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

// modify values in stiffness matrix based on Dirichlet Bcs
    for (unsigned d=0; d<NDIM; ++d)
    {
      for (unsigned side=0; side<2; ++side)
      {
        if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC)
        { // we do not need to do anything for Neumann BCs

// fetch number of nodes on face of element
          unsigned nsl = side==0 ?
            nodalBasis->getNumSurfLowerNodes(d) : nodalBasis->getNumSurfUpperNodes(d);

// allocate space for mapping
          std::vector<int> lgSurfMap(nsl);

          double dv = bc[d][side].value;
// create region to loop over side
          Lucee::Region<NDIM, int> defRgnG = side==0 ?
            globalRgn.resetBounds(d, globalRgn.getLower(d), globalRgn.getLower(d)+1) :
            globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// only update if we are on the correct ranks
          Lucee::Region<NDIM, int> defRgn = defRgnG.intersect(localRgn);

// loop, modifying stiffness matrix
          Lucee::RowMajorSequencer<NDIM> seq(defRgn);
          while (seq.step())
          {
            seq.fillWithIndex(idx);
// set index into element basis
            nodalBasis->setIndex(idx);
// get surface nodes -> global mapping
            if (side == 0)
              nodalBasis->getSurfLowerLocalToGlobal(d, lgSurfMap);
            else
              nodalBasis->getSurfUpperLocalToGlobal(d, lgSurfMap);

// reset corresponding rows (Note that some rows may be reset more
// than once. This should not be a problem, though might make the
// setup phase a bit slower).
            MatZeroRows(stiffMat, nsl, &lgSurfMap[0], 1.0);

// now insert row numbers with corresponding values into map for use
// in the update method
            for (unsigned r=0; r<nsl; ++r)
              rowBcValues[lgSurfMap[r]] = dv;
          }
        }
      }
    }

// Begin process of modification to handle periodic BCs. The code in
// the following two loops basically does the following. It modifies
// the stiffness matrix so that the nodes on the upper boundary in
// each periodic direction is identified with the corresponding node
// on the lower boundary. Further, the stiffness matrix rows
// corresponding to the nodes on the lower boundary in each periodic
// direction is also modified to take into account the contribution
// from the layer of nodes just inside the corresponding upper
// boundary.
    for (unsigned d=0; d<NDIM; ++d)
    {
      if (periodicFlgs[d] == true)
      {
// fetch number of nodes on face of element
        unsigned nsl =  nodalBasis->getNumSurfUpperNodes(d);

// space for mappings
        std::vector<int> lgUpperSurfMap(nsl), lgLowerSurfMap(nsl);
// space for local node numbers on faces
        std::vector<int> lgLocalNodeNum(nsl), lgMapMod(nlocal);
// space for stiffness mods
        std::vector<int> stiffModRowIdx(nsl);
        std::vector<double> modVals(nlocal*nlocal);

// create region to loop over side
        Lucee::Region<NDIM, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// only update if we are on the correct ranks
        Lucee::Region<NDIM, int> defRgn = defRgnG.intersect(localRgn);

// loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<NDIM> seq(defRgn);
        while (seq.step())
        {
          seq.fillWithIndex(idx);

// set index into element basis
          nodalBasis->setIndex(idx);

// get stiffness matrix for modification of the terms on the lower boundary
          nodalBasis->getStiffnessMatrix(localStiff);
// construct arrays for passing into Petsc
          for (unsigned k=0; k<nlocal; ++k)
          {
            for (unsigned m=0; m<nlocal; ++m)
              vals[nlocal*k+m] = -localStiff(k,m); // Default PetSc layout is row-major
          }

// Now compute correct locations in global stiffness matrix to add
// these values. This code looks very bizarre as the logic for
// modifying the stiffness matrix is not trivially. Basically, one
// needs to account for the wrapped node identification for the rows
// living on the lower edges while not touching the rows corresponding
// the rows on the upper edges.

// get local node number of upper face (this is on upper edge)
          nodalBasis->getSurfUpperNodeNums(d, lgLocalNodeNum);
// get local -> global mapping (this is on upper edge)
          nodalBasis->getLocalToGlobal(lgMap);
          nodalBasis->getLocalToGlobal(lgMapMod); // yes, we get this twice

// reset index to point to corresponding cell on lower edge
          idx[d] = 0;
// set index into element basis
          nodalBasis->setIndex(idx);
// get local node number of lower face (this is on lower edge)
          nodalBasis->getSurfLowerLocalToGlobal(d, lgLowerSurfMap);

// modify appropriate entries and copy over non-zero contributions
          for (unsigned k=0; k<nsl; ++k)
            lgMapMod[lgLocalNodeNum[k]-1] = lgLowerSurfMap[k];

// zero out contribution
          for (unsigned k=0; k<nlocal; ++k) modVals[k] = 0.0;
// only make contributions to those rows which correspond to those
// nodes on the lower edge
          for (unsigned k=0; k<nlocal; ++k)
          {
            for (unsigned m=0; m<nlocal; ++m)
            {
              if (lgMapMod[k] == lgMap[k])
                modVals[nlocal*k+m] = 0.0;
              else
                modVals[nlocal*k+m] = vals[nlocal*k+m];
            }
          }

// insert into global stiffness matrix
          MatSetValues(stiffMat, nlocal, &lgMapMod[0], nlocal, &lgMapMod[0],
            &modVals[0], ADD_VALUES);
        }
      }
    }
    
// reassemble matrix after modification
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

// NOTE: This second loop is needed even though it is essentially the
// same as the previous one as Petsc does not allow to call
// MatZeroRows and MatSetValues without an intervening calls to
// MatAssemblyBegin/MatAssemblyEnd. So if Petsc was not so annoying
// this extra mess would not be needed. (Ammar, April 4 2012).
    for (unsigned d=0; d<NDIM; ++d)
    {
      if (periodicFlgs[d] == true)
      {
// fetch number of nodes on face of element
        unsigned nsl =  nodalBasis->getNumSurfUpperNodes(d);

// space for mappings
        std::vector<int> lgSurfMap(nsl), lgLowerSurfMap(nsl);
// space for local node numbers on faces
        std::vector<int> lgLocalNodeNum(nsl), lgLowerLocalNodeNum(nsl);

// create region to loop over side
        Lucee::Region<NDIM, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// only update if we are on the correct ranks
        Lucee::Region<NDIM, int> defRgn = defRgnG.intersect(localRgn);

// loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<NDIM> seq(defRgn);
        while (seq.step())
        {
          seq.fillWithIndex(idx);

// set index into element basis
          nodalBasis->setIndex(idx);
// get surface nodes -> global mapping
          nodalBasis->getSurfUpperLocalToGlobal(d, lgSurfMap);

// reset corresponding rows (Note that some rows may be reset more
// than once. This should not be a problem, though might make the
// setup phase a bit slower).
          MatZeroRows(stiffMat, nsl, &lgSurfMap[0], 0.0);

// now insert row numbers with a 0.0 as corresponding source to ensure
// this point is identified with its periodic image on the lower
// boundary.
          for (unsigned r=0; r<nsl; ++r)
            rowBcValues[lgSurfMap[r]] = 0.0;

// compute corresponding node numbers on lower edge
          idx[d] = 0;
// set index into element basis
          nodalBasis->setIndex(idx);
// get surface nodes -> global mapping
          nodalBasis->getSurfLowerLocalToGlobal(d, lgLowerSurfMap);

// insert node numbers into map
          for (unsigned r=0; r<nsl; ++r)
            periodicNodeMap[lgSurfMap[r]] = lgLowerSurfMap[r];
        }
      }
    }
    
// reassemble matrix after modification
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

// list of values to force periodicity
    double periodicVals[2] = {-1, 1};
    int resetRow[1]; // row to reset
    int resetCol[2]; // col to reset

// make sweep to ensure periodicity by identifying upper edge nodes
// with lower edge nodes (periodicNodeMap holds this data)
    std::map<int, int>::const_iterator rItr
      = periodicNodeMap.begin();
    for ( ; rItr != periodicNodeMap.end(); ++rItr)
    {
      resetRow[0] = rItr->first;
      resetCol[0] = rItr->first;
      resetCol[1] = rItr->second;
      MatSetValues(stiffMat, 1, resetRow, 2, resetCol,
        periodicVals, INSERT_VALUES);
    }

// reassemble matrix after modification
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

// if all directions are periodic, set bottom left to 0.0 to avoid
// singular matrix
    if (allPeriodic)
    {
// lower-left
      int zeroRow[1] = {0};
      MatZeroRows(stiffMat, 1, zeroRow, 1.0);
// also zero out the source
      rowBcValues[0] = 0.0;
    }

// reassemble matrix after modification
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

//     PetscViewer lab;
//     PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix", &lab);
//     PetscViewerSetFormat(lab, PETSC_VIEWER_ASCII_DENSE);
//     MatView(stiffMat, lab);

//  finalize assembly
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

// create duplicate to store initial guess
    VecDuplicate(globalSrc, &initGuess);

    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, stiffMat, stiffMat, DIFFERENT_NONZERO_PATTERN);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetFromOptions(ksp);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FemPoissonStructUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input/output fields
    const Lucee::Field<NDIM, double>& src = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& sol = this->getOut<Lucee::Field<NDIM, double> >(0);

// global region for grid
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
// create sequencer for looping over local box
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// local mass matrix
    Lucee::Matrix<double> localMass(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// storage for computing source contribution
    std::vector<double> localSrc(nlocal), localMassSrc(nlocal);

// pointers
    Lucee::ConstFieldPtr<double> srcPtr = src.createConstPtr();

// create sequencer for looping over local box
    Lucee::RowMajorSequencer<NDIM> seq(grid.getLocalRegion());
    int idx[NDIM];

    double vol = grid.getComputationalSpace().getVolume();
    double intSrcVol = 0.0;
// if both directions are periodic, we need to adjust source to ensure
// solvability of the equations
    if (allPeriodic)
      intSrcVol = getFieldIntegral(src, srcNodesShared)/vol;

// clear out existing stuff in source vector: this is required
// otherwise successive calls to advance() will accumulate into source
// from prevous calls, which is of course not what we want.
    VecSet(globalSrc, 0.0);

// loop, creating RHS (source terms)
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// set index into element basis
      nodalBasis->setIndex(idx);

// get local mass matrix
      nodalBasis->getMassMatrix(localMass);
// get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

      if (srcNodesShared)
      {
// extract source at each node from field
        nodalBasis->extractFromField(src, localSrc);
      }
      else
      {
        src.setPtr(srcPtr, idx);
// if nodes are not shared simply copy over data
        for (unsigned k=0; k<nlocal; ++k)
          localSrc[k] = srcPtr[k];
      }

// adjust for periodic BCs in all directions
      for (unsigned k=0; k<nlocal; ++k)
        localSrc[k] += -intSrcVol;

// evaluate local mass matrix times local source
      for (unsigned k=0; k<nlocal; ++k)
      {
        localMassSrc[k] = 0.0;
        for (unsigned m=0; m<nlocal; ++m)
          localMassSrc[k] += localMass(k,m)*localSrc[m];
      }
// accumulate it into Petsc vector
      VecSetValues(globalSrc, nlocal, &lgMap[0], &localMassSrc[0], ADD_VALUES);
    }

// finish assembly of RHS
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

// Now loop over each periodic direction, adding contribution from
// last layer of cells on upper edges to lower edge nodes
    for (unsigned d=0; d<NDIM; ++d)
    {
      if (periodicFlgs[d] == true)
      {
// fetch number of nodes on face of element
        unsigned nsl =  nodalBasis->getNumSurfUpperNodes(d);

// space for mappings
        std::vector<int> lgLowerSurfMap(nsl);
// space for local node numbers on faces
        std::vector<int> lgLocalNodeNum(nsl);
// space for passing into Petsc
        std::vector<double> localMassMod(nsl);

// create region to loop over side
        Lucee::Region<NDIM, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// only update if we are on the correct ranks
        Lucee::Region<NDIM, int> defRgn = defRgnG.intersect(localRgn);

// loop, modifying source vector
        Lucee::RowMajorSequencer<NDIM> seq(defRgn);
        while (seq.step())
        {
          seq.fillWithIndex(idx);

// set index into element basis
          nodalBasis->setIndex(idx);
// get local node number of upper face
          nodalBasis->getSurfUpperNodeNums(d, lgLocalNodeNum);

// get local mass matrix
          nodalBasis->getMassMatrix(localMass);

          if (srcNodesShared)
          {
// extract source at each node from field
            nodalBasis->extractFromField(src, localSrc);
          }
          else
          {
            src.setPtr(srcPtr, idx);
// if nodes are not shared simply copy over data
            for (unsigned k=0; k<nlocal; ++k)
              localSrc[k] = srcPtr[k];
          }

// adjust for periodic BCs in all directions
          for (unsigned k=0; k<nlocal; ++k)
            localSrc[k] += -intSrcVol;

// evaluate local mass matrix times local source
          for (unsigned k=0; k<nlocal; ++k)
          {
            localMassSrc[k] = 0.0;
            for (unsigned m=0; m<nlocal; ++m)
              localMassSrc[k] += localMass(k,m)*localSrc[m];
          }

// just copy appropriate data over for passing into PetSc
          for (unsigned k=0; k<nsl; ++k)
            localMassMod[k] = localMassSrc[lgLocalNodeNum[k]-1];

// compute corresponding node numbers on lower edge
          idx[d] = 0;
// set index into element basis
          nodalBasis->setIndex(idx);
// get surface nodes -> global mapping
          nodalBasis->getSurfLowerLocalToGlobal(d, lgLowerSurfMap);

// accumulate it into Petsc vector
          VecSetValues(globalSrc, nsl, &lgLowerSurfMap[0], &localMassMod[0], ADD_VALUES);
        }
      }
    }

// finish assembly of RHS
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

    int resetRow[1];
    double resetVal[1];
// reset source to apply boundary conditions
    std::map<int, double>::const_iterator rItr
      = rowBcValues.begin();
    for ( ; rItr != rowBcValues.end(); ++rItr)
    {
      resetRow[0] = rItr->first; // row index
      resetVal[0] = rItr->second; // value
      VecSetValues(globalSrc, 1, resetRow, resetVal, INSERT_VALUES);
    }

// reassemble RHS after application of Dirichlet Bcs
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

//     PetscViewer lab;
//     PetscViewerASCIIOpen(PETSC_COMM_WORLD, "vector", &lab);
//     PetscViewerSetFormat(lab, PETSC_VIEWER_DEFAULT);
//     VecView(globalSrc, lab);

// copy solution for use as initial guess in KSP solve
    PetscScalar *ptGuess;
    unsigned count = 0;
    VecGetArray(initGuess, &ptGuess);
    nodalBasis->copyAllDataFromField(sol, ptGuess);
    VecRestoreArray(initGuess, &ptGuess);

// now solve linear system (initGuess will contain solution)
    KSPSolve(ksp, globalSrc, initGuess);
// check if solver converged
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    int itNum;
    KSPGetIterationNumber(ksp, &itNum);

// construct message to send back to Lua
    std::ostringstream msgStrm;
    bool status = true;
    if (reason < 0) 
    {
      status = false;
      msgStrm << FemPoissonStructUpdater<NDIM>::id << ": KSPSolve failed!";
      msgStrm << " Petsc reason code was " << reason << ".";
    }
    else
    {
      msgStrm << FemPoissonStructUpdater<NDIM>::id << ": KSPSolve converged.";
    }
    msgStrm << " Number of iterations " << itNum;

// copy solution from PetSc array to solution field
    PetscScalar *ptSol;
    VecGetArray(initGuess, &ptSol);
    nodalBasis->copyAllDataToField(ptSol, sol);
    VecRestoreArray(initGuess, &ptSol);

    return Lucee::UpdaterStatus(status, std::numeric_limits<double>::max(),
      msgStrm.str());
  }

  template <unsigned NDIM>
  void
  FemPoissonStructUpdater<NDIM>::declareTypes()
  {
// takes one input (source terms)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// returns one output, solution
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  typename FemPoissonStructUpdater<NDIM>::FemPoissonBcData
  FemPoissonStructUpdater<NDIM>::getBcData(const Lucee::LuaTable& bct) const
  {
    FemPoissonBcData bcData;
    if (bct.getString("T") == "D")
      bcData.type = DIRICHLET_BC;
    else if ((bct.getString("T") == "N"))
      bcData.type = NEUMANN_BC;
    else
    {
      Lucee::Except lce(
        "FemPoissonStructUpdater::readInput: Must specify one of \"D\" or \"N\". ");
      lce << "Specified \"" << bct.getString("T") << " instead";
      throw lce;
    }
    bcData.value = bct.getNumber("V");
    bcData.isSet = true;
    
    return bcData;
  }

  template <unsigned NDIM>
  double
  FemPoissonStructUpdater<NDIM>::getFieldIntegral(const Lucee::Field<NDIM, double>& fld, 
    bool shareFlag)
  {
    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<double> weights(nlocal), localFld(nlocal);
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// create pointers
    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();

    Lucee::RowMajorSequencer<NDIM> seq(grid.getLocalRegion());
    int idx[NDIM];
    double fldInt = 0.0;
// loop, accumulating integral in each cell
    while (seq.step())
    {
      seq.fillWithIndex(idx);

// set index into element basis
      nodalBasis->setIndex(idx);

// get quadrature weights
      nodalBasis->getWeights(weights);

      if (shareFlag)
      {
// extract source at each node from field
        nodalBasis->extractFromField(fld, localFld);
      }
      else
      {
        fld.setPtr(fldPtr, idx);
// if nodes are not shared simply copy over data
        for (unsigned k=0; k<nlocal; ++k)
          localFld[k] = fldPtr[k];
      }

// compute contribition from this cell
      for (unsigned k=0; k<nlocal; ++k)
      {
        fldInt += weights[k]*localFld[k];      
      }
    }

    double netFldInt = fldInt;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
// sum across all processors
    comm->allreduce(1, &fldInt, &netFldInt, TX_SUM);

    return netFldInt;
  }

// instantiations
  template class FemPoissonStructUpdater<1>;
  template class FemPoissonStructUpdater<2>;
  template class FemPoissonStructUpdater<3>;
}

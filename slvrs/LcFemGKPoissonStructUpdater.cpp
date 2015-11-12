/**
 * @file	LcFemGKPoissonStructUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with FEM scheme on a structured grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFemGKPoissonStructUpdater.h>
#include <LcGlobals.h>
#include <LcMatrix.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

// txbase includes
#include <TxCommBase.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <set>
#include <sstream>
#include <vector>
#include <limits>

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

  template <> const char *FemGKPoissonStructUpdater<2>::id = "FemGKPoisson2D";
  template <> const char *FemGKPoissonStructUpdater<3>::id = "FemGKPoisson3D";

  template <unsigned NDIM>
  FemGKPoissonStructUpdater<NDIM>::FemGKPoissonStructUpdater()
    : Lucee::UpdaterIfc()
  {
// this flag is needed to ensure PetSc arrays are only deleted if
// update() method is called at least once.
    runOnce = false;
  }

  template <unsigned NDIM>
  FemGKPoissonStructUpdater<NDIM>::~FemGKPoissonStructUpdater()
  {
    if (runOnce)
    {
      MatDestroy(stiffMat);
      VecDestroy(globalSrc);
      VecDestroy(initGuess);
      //VecDestroy(localData);
      //ISDestroy(is);
    }
  }

  template <unsigned NDIM>
  void
  FemGKPoissonStructUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except(
        "FemGKPoissonStructUpdater::readInput: Must specify element to use using 'basis'");

// check if source nodes are shared
    srcNodesShared = true;
    if (tbl.hasBool("sourceNodesShared"))
      srcNodesShared = tbl.getBool("sourceNodesShared");

// check if solution nodes are shared
    solNodesShared = true;
    if (tbl.hasBool("solutionNodesShared"))
      solNodesShared = tbl.getBool("solutionNodesShared");

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
          throw Lucee::Except("FemGKPoissonStructUpdater:readInput: Incorrect 'periodicDirs'");
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
          Lucee::Except lce("FemGKPoissonStructUpdater:: Must specify BCs on each side");
          throw lce;
        }
      }
      else
      {
// ensure no BCs are specified
        if (bc[dir][0].isSet || bc[dir][1].isSet)
        {
          Lucee::Except lce("FemGKPoissonStructUpdater:: Cannot specify BCs");
          lce << " if a direction is periodic";
          throw lce;
        }
      }
    }

    modifierConstant = 0.0;
// read in modifier constant (if any)
    if (tbl.hasNumber("modifierConstant"))
      modifierConstant = tbl.getNumber("modifierConstant");

// Adjust source if doing periodic BCs and no modifier constant specified
    adjustSource = false;
    if (allPeriodic && (modifierConstant==0) )
      adjustSource = true;

// check if stiffness matrix should be written out
    writeMatrix = false;
    if (tbl.hasBool("writeStiffnessMatrix"))
      writeMatrix = tbl.getBool("writeStiffnessMatrix");
  }

  template <unsigned NDIM>
  void
  FemGKPoissonStructUpdater<NDIM>::initialize()
  {
//#define DMSG(s) std::cout << "Rank " << this->getComm()->getRank() << ": " << s << std::endl;
#define DMSG(s) ;

    DMSG("Inside initialize");

    Lucee::UpdaterIfc::initialize();

    unsigned nglobal = nodalBasis->getNumGlobalNodes();
    unsigned nlocal = nodalBasis->getNumNodes();
// map to hold periodicially identified nodes
    std::map<int, int> periodicNodeMap;
#ifdef HAVE_MPI
    TxMpiBase *comm = static_cast<TxMpiBase*>(this->getComm());
#endif

#ifdef HAVE_MPI
    int nz = nodalBasis->getNumNodes()*(std::pow(2.0, 1.0*NDIM)+1);
    MatCreateMPIAIJ(comm->getMpiComm(), PETSC_DECIDE, PETSC_DECIDE, nglobal, nglobal, 
      nz, PETSC_NULL, 
      nz, PETSC_NULL,
      &stiffMat);
#else
// Explicit initialization of stiffness matrix speeds up
// initialization tremendously.
    int nz = nodalBasis->getNumNodes()*(std::pow(2.0, 1.0*NDIM)+1);
    MatCreateSeqAIJ(PETSC_COMM_SELF, nglobal, nglobal, nz, PETSC_NULL, &stiffMat);
#endif
    MatSetFromOptions(stiffMat);

// now create vector to store source: MatGetVecs ensures that the
// parallel layout is the same as the stiffness matrix.
    MatGetVecs(stiffMat, &globalSrc, PETSC_NULL);
    VecSetFromOptions(globalSrc);

    DMSG("Matrix and vectors allocated");

    Lucee::Matrix<double> localStiff(nlocal, nlocal), localMass(nlocal, nlocal);
    std::vector<int> lgMap(nlocal);
    std::vector<PetscScalar> vals(nlocal*nlocal);

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion(); 
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    int idx[NDIM];

// loop, creating global stiffness matrix
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);

      nodalBasis->getPerpStiffnessMatrix(localStiff);
      nodalBasis->getMassMatrix(localMass);
// construct arrays for passing into Petsc
      for (unsigned k=0; k<nlocal; ++k)
      {
        for (unsigned m=0; m<nlocal; ++m)
          vals[nlocal*k+m] = -localStiff(k,m)+modifierConstant*localMass(k,m); //  // Default PetSc layout is row-major
      }

      nodalBasis->getLocalToGlobal(lgMap);
// insert into global stiffness matrix, adding them to existing value
      MatSetValues(stiffMat, nlocal, &lgMap[0], nlocal, &lgMap[0],
        &vals[0], ADD_VALUES);
    }

    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

    DMSG("Basic stiffness matrix computed");

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

// loop, modifying stiffness matrix NOTE: We need are running this
// loop on all processors, even those that don't own the boundary
// nodes. The reason is that PetSc expects this.
          Lucee::RowMajorSequencer<NDIM> seq(defRgnG);
          while (seq.step())
          {
            seq.fillWithIndex(idx);
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

    DMSG("Modification for Dirichlet BCs completed");

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

        std::vector<int> lgUpperSurfMap(nsl), lgLowerSurfMap(nsl);
        std::vector<int> lgLocalNodeNum(nsl), lgMapMod(nlocal);
        std::vector<int> stiffModRowIdx(nsl);
        std::vector<double> modVals(nlocal*nlocal);

// create region to loop over side
        Lucee::Region<NDIM, int> defRgn = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// region for ensuring that we only update the stiffness matrix on
// appropriate rank
        Lucee::Region<NDIM, int> defRgnUp = defRgn.intersect(localRgn);

// loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<NDIM> seq(defRgn);
        while (seq.step())
        {
          seq.fillWithIndex(idx);
          nodalBasis->setIndex(idx);

// this flag is needed to ensure we don't update the stiffness matrix twice
          bool isIdxLocal = defRgnUp.isInside(idx);

          nodalBasis->getPerpStiffnessMatrix(localStiff);
          nodalBasis->getMassMatrix(localMass);          
// construct arrays for passing into Petsc
          for (unsigned k=0; k<nlocal; ++k)
          {
            for (unsigned m=0; m<nlocal; ++m)
              vals[nlocal*k+m] = -localStiff(k,m)+modifierConstant*localMass(k,m); // Default PetSc layout is row-major
          }

// Now compute correct locations in global stiffness matrix to add
// these values. This code looks very bizarre as the logic for
// modifying the stiffness matrix is not trivial. Basically, one needs
// to account for the wrapped node identification for the rows living
// on the lower edges while not touching the rows corresponding the
// rows on the upper edges.

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
            lgMapMod[lgLocalNodeNum[k]] = lgLowerSurfMap[k];

// zero out contribution
          //for (unsigned k=0; k<nlocal; ++k) modVals[k] = 0.0;
          for (unsigned k=0; k<nlocal*nlocal; ++k) modVals[k] = 0.0;

// the following check is needed to ensure that the periodic BC
// modifications to stiffness matrix are not applied twice
          if (isIdxLocal)
          {
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

    DMSG("Phase I of periodic BCs completed");

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

        std::vector<int> lgSurfMap(nsl), lgLowerSurfMap(nsl);
        std::vector<int> lgLocalNodeNum(nsl), lgLowerLocalNodeNum(nsl);

// create region to loop over side
        Lucee::Region<NDIM, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// only update if we are on the correct ranks
        Lucee::Region<NDIM, int> defRgn = defRgnG;//.intersect(localRgn);

// loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<NDIM> seq(defRgn);
        while (seq.step())
        {
          seq.fillWithIndex(idx);
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

    DMSG("Phase II of periodic BCs completed");

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
    if (adjustSource) //if (allPeriodic)
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

    if (writeMatrix)
    {
      std::string outName = Loki::SingletonHolder<Lucee::Globals>
        ::Instance().outPrefix + "-poisson-stiffnessMatrix";
      PetscViewer lab;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, outName.c_str(), &lab);
      PetscViewerSetFormat(lab, PETSC_VIEWER_ASCII_MATLAB);
      MatView(stiffMat, lab);
    }

//  finalize assembly
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

// create duplicate to store initial guess (this will also contain
// solution)
    VecDuplicate(globalSrc, &initGuess);

// Create index set to copy data from parallel Petsc vector to locally
// required data. 
//
// NOTE: This is required to ensure that the solution is available
// with same parallel layout as expected by Gkeyll. Otherwise, things
// go hay-wire as Gkeyll and Petsc use slightly different parallel
// layouts, causing mayhem and chaos.
    std::vector<PetscInt> vecIs;
    Lucee::RowMajorSequencer<NDIM> seqIs(localRgn);
    while (seqIs.step())
    {
      seqIs.fillWithIndex(idx);
      nodalBasis->setIndex(idx);
      nodalBasis->getLocalToGlobal(lgMap);
      for (unsigned k=0; k<nlocal; ++k)
        vecIs.push_back(lgMap[k]);
    }

    PetscInt numIs = vecIs.size();
#ifdef HAVE_MPI
    /*Petsc 3.2+ ISCreateGeneral(comm->getMpiComm(), numIs, &vecIs[0], PETSC_COPY_VALUES, &is);*/
    ISCreateGeneral(comm->getMpiComm(), numIs, &vecIs[0], &is);
    VecCreateMPI(comm->getMpiComm(), numIs, PETSC_DETERMINE, &localData);
#else
    /*Petsc 3.2+ ISCreateGeneral(PETSC_COMM_SELF, numIs, &vecIs[0], PETSC_COPY_VALUES, &is);*/
    ISCreateGeneral(PETSC_COMM_SELF, numIs, &vecIs[0], &is);
    VecCreateSeq(PETSC_COMM_SELF, numIs, &localData);
#endif

// create a scatter context to get data onto local processors
    VecScatterCreate(initGuess, is, localData, PETSC_NULL, &vecSctr);

    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, stiffMat, stiffMat, DIFFERENT_NONZERO_PATTERN);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetFromOptions(ksp);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FemGKPoissonStructUpdater<NDIM>::update(double t)
  {
// set flag to indicate we should delete PetSc vectors/matrices
    runOnce = true;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& src = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& sol = this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    unsigned nlocal = nodalBasis->getNumNodes();

    Lucee::Matrix<double> localMass(nlocal, nlocal);
    std::vector<int> lgMap(nlocal);
    std::vector<double> localSrc(nlocal), localMassSrc(nlocal);
    Lucee::ConstFieldPtr<double> srcPtr = src.createConstPtr();

    Lucee::RowMajorSequencer<NDIM> seq(grid.getLocalRegion());
    int idx[NDIM];

    double vol = grid.getComputationalSpace().getVolume();
    double intSrcVol = 0.0;
// if all directions are periodic, we need to adjust source to ensure
// solvability of the equations
    if (adjustSource) //if (allPeriodic)
      intSrcVol = getFieldIntegral(src, srcNodesShared)/vol;

// clear out existing stuff in source vector: this is required
// otherwise successive calls to advance() will accumulate into source
// from prevous calls, which is of course not what we want.
    VecSet(globalSrc, 0.0);

// loop, creating RHS (source terms)
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);

      nodalBasis->getMassMatrix(localMass);
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

        std::vector<int> lgLowerSurfMap(nsl);
        std::vector<int> lgLocalNodeNum(nsl);
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
          nodalBasis->setIndex(idx);
          nodalBasis->getSurfUpperNodeNums(d, lgLocalNodeNum);
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
            localMassMod[k] = localMassSrc[lgLocalNodeNum[k]];

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

    if (writeMatrix)
    {
      std::string outName = Loki::SingletonHolder<Lucee::Globals>
        ::Instance().outPrefix + "-poisson-globalSrc";
      PetscViewer lab;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, outName.c_str(), &lab);
      PetscViewerSetFormat(lab, PETSC_VIEWER_ASCII_MATLAB);
      VecView(globalSrc, lab);
    }

// copy solution for use as initial guess in KSP solve
    copyFromGkeyllField(sol, initGuess);

// now solve linear system (initGuess will contain solution)
    KSPSolve(ksp, globalSrc, initGuess);
// check if solver converged
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    int itNum;
    KSPGetIterationNumber(ksp, &itNum);
    double resNorm;
    KSPGetResidualNorm(ksp, &resNorm);

    if (writeMatrix)
    {
      std::string outName = Loki::SingletonHolder<Lucee::Globals>
        ::Instance().outPrefix + "-poisson-solution";
      PetscViewer lab;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, outName.c_str(), &lab);
      PetscViewerSetFormat(lab, PETSC_VIEWER_ASCII_MATLAB);
      VecView(initGuess, lab);
    }

// construct message to send back to Lua
    std::ostringstream msgStrm;
    bool status = true;
    if (reason < 0) 
    {
      status = false;
      msgStrm << FemGKPoissonStructUpdater<NDIM>::id << ": KSPSolve failed!";
      msgStrm << " Petsc reason code was " << reason << ".";
    }
    else
    {
      msgStrm << FemGKPoissonStructUpdater<NDIM>::id << ": KSPSolve converged.";
    }
    msgStrm << " Number of iterations " << itNum
            << ". Final residual norm was " << resNorm;
    //std::cout << msgStrm.str() << std::endl;

// copy solution from PetSc array to solution field
    copyFromPetscField(initGuess, sol);

    return Lucee::UpdaterStatus(status, std::numeric_limits<double>::max(),
      msgStrm.str());
  }

  template <unsigned NDIM>
  void
  FemGKPoissonStructUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  typename FemGKPoissonStructUpdater<NDIM>::FemPoissonBcData
  FemGKPoissonStructUpdater<NDIM>::getBcData(const Lucee::LuaTable& bct) const
  {
    FemPoissonBcData bcData;
    if (bct.getString("T") == "D")
      bcData.type = DIRICHLET_BC;
    else if ((bct.getString("T") == "N"))
      bcData.type = NEUMANN_BC;
    else
    {
      Lucee::Except lce(
        "FemGKPoissonStructUpdater::readInput: Must specify one of \"D\" or \"N\". ");
      lce << "Specified \"" << bct.getString("T") << " instead";
      throw lce;
    }
    bcData.value = bct.getNumber("V");
    bcData.isSet = true;
    
    return bcData;
  }

  template <unsigned NDIM>
  double
  FemGKPoissonStructUpdater<NDIM>::getFieldIntegral(const Lucee::Field<NDIM, double>& fld, 
    bool shareFlag)
  {
    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<double> weights(nlocal), localFld(nlocal);
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    Lucee::RowMajorSequencer<NDIM> seq(grid.getLocalRegion());
    int idx[NDIM];
    double fldInt = 0.0;
// loop, accumulating integral in each cell
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);
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

      for (unsigned k=0; k<nlocal; ++k)
        fldInt += weights[k]*localFld[k];      
    }

    double netFldInt = fldInt;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &fldInt, &netFldInt, TX_SUM);

    return netFldInt;
  }

  template <unsigned NDIM>
  void
  FemGKPoissonStructUpdater<NDIM>::copyFromPetscField(Vec ptFld, Lucee::Field<NDIM, double>& gkFld)
  {
    PetscScalar *ptFldPtr;
    if (solNodesShared)
    {
      VecGetArray(ptFld, &ptFldPtr);
      nodalBasis->copyAllDataToField(ptFldPtr, gkFld);
      VecRestoreArray(ptFld, &ptFldPtr);
    }
    else
    {
// get data local to this processor (this is needed as Petsc parallel
// layout is slightly different than Gkeyll layout. Most data is
// already local and does not require communication.)
      VecScatterBegin(vecSctr, ptFld, localData, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(vecSctr, ptFld, localData, INSERT_VALUES, SCATTER_FORWARD);

// copy data over to Gkeyll field
      const Lucee::StructuredGridBase<NDIM>& grid 
        = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

      unsigned count = 0;
      unsigned nlocal = nodalBasis->getNumNodes();
      Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
      Lucee::FieldPtr<double> gkPtr = gkFld.createPtr();
      Lucee::RowMajorSequencer<NDIM> seq(localRgn);
      int idx[NDIM];

      VecGetArray(localData, &ptFldPtr);
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        gkFld.setPtr(gkPtr, idx);
        for (unsigned k=0; k<nlocal; ++k)
          gkPtr[k] = ptFldPtr[count++];
      }
      VecRestoreArray(localData, &ptFldPtr);
    }
  }

  template <unsigned NDIM>
  void
  FemGKPoissonStructUpdater<NDIM>::copyFromGkeyllField(const Lucee::Field<NDIM, double>& gkFld, Vec ptFld)
  {
    PetscScalar *ptFldPtr;
    if (solNodesShared)
    {
      VecGetArray(ptFld, &ptFldPtr);
      nodalBasis->copyAllDataFromField(gkFld, ptFldPtr);
      VecRestoreArray(ptFld, &ptFldPtr);
    }
    else
    {
// copy data from Gkeyll field
      const Lucee::StructuredGridBase<NDIM>& grid 
        = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

      unsigned count = 0;
      unsigned nlocal = nodalBasis->getNumNodes();
      Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
      Lucee::ConstFieldPtr<double> gkPtr = gkFld.createConstPtr();
      Lucee::RowMajorSequencer<NDIM> seq(localRgn);
      int idx[NDIM];

      VecGetArray(localData, &ptFldPtr);
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        gkFld.setPtr(gkPtr, idx);
        for (unsigned k=0; k<nlocal; ++k)
          ptFldPtr[count++] = gkPtr[k];
      }
      VecRestoreArray(localData, &ptFldPtr);

// get data local to this processor (this is needed as Petsc parallel
// layout is slightly different than Gkeyll layout. Most data is
// already local and does not require communication.)
      VecScatterBegin(vecSctr, localData, ptFld, INSERT_VALUES, SCATTER_REVERSE);
      VecScatterEnd(vecSctr, localData, ptFld, INSERT_VALUES, SCATTER_REVERSE);
    }
  }

// instantiations
  template class FemGKPoissonStructUpdater<2>;
  template class FemGKPoissonStructUpdater<3>;
}

// The only place solution appears is when setting an initial
// condition for iterative solver and when copying solution back to
// Gkeyll fields. This is now the remaining problem in getting the
// correct solution in parallel.
//
// One thing I don't understand is how the local->global mapping can
// work.
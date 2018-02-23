/**
 * @file	LcElectrostaticContPhiUpdater.cpp
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcElectrostaticContPhiUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// from the contfromdiscontupdater
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
#include <limits>

// for cutoff velocities
#include <LcDynVector.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *ElectrostaticContPhiUpdater::id = "ElectrostaticContPhiUpdater";
  const unsigned NDIM = 1;

  ElectrostaticContPhiUpdater::ElectrostaticContPhiUpdater()
    : UpdaterIfc()
  {
  }

  ElectrostaticContPhiUpdater::~ElectrostaticContPhiUpdater()
  {
#if PETSC_VERSION_GE(3,6,0)
    MatDestroy(&stiffMat);
    VecDestroy(&globalSrc);
    VecDestroy(&initGuess);
#else
    MatDestroy(stiffMat);
    VecDestroy(globalSrc);
    VecDestroy(initGuess);
#endif
  }

  void 
  ElectrostaticContPhiUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("ElectrostaticContPhiUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerpTimesRho"))
      kPerpTimesRho = tbl.getNumber("kPerpTimesRho");
    else
      throw Lucee::Except("ElectrostaticContPhiUpdater::readInput: Must specify kPerpTimesRho");

    if (tbl.hasNumber("Te0"))
      Te0 = tbl.getNumber("Te0");
    else
      throw Lucee::Except("ElectrostaticContPhiUpdater::readInput: Must specify Te0");

    useCutoffVelocities = false;
    if (tbl.hasBool("useCutoffVelocities"))
      useCutoffVelocities = tbl.getBool("useCutoffVelocities");

    // force all directions to be periodic
    for (int i = 0; i < NDIM; i++) 
      periodicFlgs[i] = true;
  }

  void 
  ElectrostaticContPhiUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // number of global nodes
    int nglobal = nodalBasis->getNumGlobalNodes();
    // number of local nodes
    int nlocal = nodalBasis->getNumNodes();
#ifdef HAVE_MPI
    TxMpiBase *comm = static_cast<TxMpiBase*>(this->getComm());
#endif
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    // local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    // global region to update
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion();

//    Lucee::RowMajorSequencer<NDIM> seq(globalRgn);
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get mass matrix and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);

    massMatrix = Eigen::MatrixXd(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);
    gaussOrdinates = Eigen::MatrixXd(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    tripleProducts = std::vector<Eigen::MatrixXd>(nlocal);
    // Create and store the triple-product basis integrals

    for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
    {
      // Initialize matrices
      tripleProducts[basisIndex]  = Eigen::MatrixXd::Zero(nlocal, nlocal);
      for (int rowIndex = 0; rowIndex < nlocal; rowIndex++)
      {
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
        {
          double integralResult = 0.0;

          for (int gaussNodeIndex = 0; gaussNodeIndex < numQuadNodes; gaussNodeIndex++)
          {
            integralResult += gaussWeights[gaussNodeIndex]*interpMatrix(gaussNodeIndex, basisIndex)*
              interpMatrix(gaussNodeIndex, rowIndex)*interpMatrix(gaussNodeIndex, colIndex);
          }

          tripleProducts[basisIndex](rowIndex, colIndex) = integralResult;
        }
      }
    }

#ifdef HAVE_MPI
    int nz = nodalBasis->getNumNodes()*(std::pow(2.0, 1.0*NDIM)+1);
#if PETSC_VERSION_GE(3,6,0)
    MatCreateAIJ(comm->getMpiComm(), PETSC_DECIDE, PETSC_DECIDE, nglobal, nglobal,
#else
    MatCreateMPIAIJ(comm->getMpiComm(), PETSC_DECIDE, PETSC_DECIDE, nglobal, nglobal,
#endif
      nz, PETSC_NULL,
      nz, PETSC_NULL,
      &stiffMat);
#else
    // Explicit initialization of stiffness matrix speeds up
    // initialization tremendously.
    //int nz = 10; // number of non-zero entries per row (WHAT SHOULD IT REALLY BE?)
    //if (nodalBasis->getNumNodes() > 4)
    //  nz = 21; // HACK FOR NOW
    int nz = nodalBasis->getNumNodes()*(std::pow(2.0, 1.0*NDIM)+1);
    MatCreateSeqAIJ(PETSC_COMM_SELF, nglobal, nglobal, nz, PETSC_NULL, &stiffMat);
#endif
    MatSetFromOptions(stiffMat);

#if PETSC_VERSION_GE(3,6,0)
// From LcFemPoissonStructUpdater.cpp: NRM 3/10/16:
// Depending on the type of BCs, we will be modifying stiffMat and adding new nonzero values.
// In later versions of PETSc, we need to set the following in order for PETSc to not throw errors.
    MatSetOption(stiffMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif
// Keep non-zero pattern when zero-ing out matrix
    MatSetOption(stiffMat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    // Finalize assembly... is this needed?
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);
    
    // create and setup vector:
    // (These work in serial)
//    VecCreate(MPI_COMM_WORLD, &globalSrc);
//    VecSetSizes(globalSrc, nglobal, PETSC_DECIDE);
//    VecSetFromOptions(globalSrc);
    // (NOTE: ACCORDING TO MIKE MCCOURT ONE SHOULD
    // USE MatGetVecs INSTEAD OF VecCreate. THIS ENSURES THAT THE PARALLEL
    // LAYOUT OF THE MATRIX & VECTOR ARE THE SAME)
#if PETSC_VERSION_GE(3,6,0)
    MatCreateVecs(stiffMat, &globalSrc, PETSC_NULL);
#else
    MatGetVecs(stiffMat, &globalSrc, PETSC_NULL);
#endif
    VecSetFromOptions(globalSrc);
    
    // finalize assembly
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);
    
    // create duplicate to store initial guess
    VecDuplicate(globalSrc, &initGuess); 

// Create index set to copy data from parallel Petsc vector to locally
// required data.
//
// NOTE: This is required to ensure that the solution is available
// with same parallel layout as expected by Gkeyll. Otherwise, things
// go hay-wire as Gkeyll and Petsc use slightly different parallel
// layouts, causing mayhem and chaos.
    std::vector<int> lgMap(nlocal);
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
# if PETSC_VERSION_GE(3,6,0)
    ISCreateGeneral(comm->getMpiComm(), numIs, &vecIs[0], PETSC_COPY_VALUES, &is);
# else
    ISCreateGeneral(comm->getMpiComm(), numIs, &vecIs[0], &is);
# endif
    VecCreateMPI(comm->getMpiComm(), numIs, PETSC_DETERMINE, &localData);

#else
# if PETSC_VERSION_GE(3,6,0)
    ISCreateGeneral(PETSC_COMM_SELF, numIs, &vecIs[0], PETSC_COPY_VALUES, &is);
# else
    ISCreateGeneral(PETSC_COMM_SELF, numIs, &vecIs[0], &is);
# endif
    VecCreateSeq(PETSC_COMM_SELF, numIs, &localData);
#endif

// create a scatter context to get data onto local processors
    VecScatterCreate(initGuess, is, localData, PETSC_NULL, &vecSctr);
    
    // Create solver
    KSPCreate(MPI_COMM_WORLD, &ksp);
#if PETSC_VERSION_GE(3,6,0)
    KSPSetOperators(ksp, stiffMat, stiffMat);
#else
    KSPSetOperators(ksp, stiffMat, stiffMat, DIFFERENT_NONZERO_PATTERN);
#endif
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetFromOptions(ksp);
  }

  Lucee::UpdaterStatus 
  ElectrostaticContPhiUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // Electron and ion densities
    const Lucee::Field<NDIM, double>& nElcIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& nIonIn = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& phiOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::ConstFieldPtr<double> nElcPtr = nElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::FieldPtr<double> phiPtr = phiOut.createPtr();

    // number of global nodes
    int nglobal = nodalBasis->getNumGlobalNodes();
    // number of local nodes
    int nlocal = nodalBasis->getNumNodes();

    // Local mass matrix
    Lucee::Matrix<double> localMassLucee(nlocal, nlocal);

    // map to hold periodicially identified nodes
    std::map<int, int> periodicNodeMap;

    // Zero out mat matrix
    MatZeroEntries(stiffMat);
    // reassemble matrix after modification
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

    // map for local indices to global indices
    std::vector<int> lgMap(nlocal);

    // storage for passing to petsc
    std::vector<PetscScalar> vals(nlocal*nlocal);

    // global region for grid
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
    // create sequencer for looping over local box
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion(); 
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    int idx[NDIM];

    // loop, creating global stiffness matrix
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // set index into element basis
      nodalBasis->setIndex(idx);
      // set index of nIon
      nIonIn.setPtr(nIonPtr, idx);

      Eigen::VectorXd nIonVec(nlocal);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        nIonVec(componentIndex) = nIonPtr[componentIndex];

      Eigen::MatrixXd phiProjections(nlocal, nlocal);

      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        Eigen::VectorXd phiProjectionsSingle = tripleProducts[basisIndex]*nIonVec;
        // Store components into rows of phiProjections matrix
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
          phiProjections(basisIndex, colIndex) = phiProjectionsSingle(colIndex);
      }

      // construct arrays for passing into Petsc
      for (int k = 0; k < nlocal; k++)
        for (int m = 0; m < nlocal; m++)
          vals[nlocal*k + m] = phiProjections(k, m); // Default PetSc layout is row-major

      // get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

      // insert into global stiffness matrix, adding them to existing value
      MatSetValues(stiffMat, nlocal, &lgMap[0], nlocal, &lgMap[0],
        &vals[0], ADD_VALUES);
    }

    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

    // Begin process of modification to handle periodic BCs. The code in
    // the following two loops basically does the following. It modifies
    // the stiffness matrix so that the nodes on the upper boundary in
    // each periodic direction is identified with the corresponding node
    // on the lower boundary. Further, the stiffness matrix rows
    // corresponding to the nodes on the lower boundary in each periodic
    // direction is also modified to take into account the contribution
    // from the layer of nodes just inside the corresponding upper
    // boundary.
    for (unsigned d = 0; d < NDIM; d++)
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
        Lucee::RowMajorSequencer<NDIM> seqSide(defRgn);
        while (seqSide.step())
        {
          seqSide.fillWithIndex(idx);

          // set index into element basis
          nodalBasis->setIndex(idx);
          // set index of nIon
          nIonIn.setPtr(nIonPtr, idx);

// this flag is needed to ensure we don't update the stiffness matrix twice
          bool isIdxLocal = defRgn.isInside(idx);

          Eigen::VectorXd nIonVec(nlocal);
          
          for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
            nIonVec(componentIndex) = nIonPtr[componentIndex];

          Eigen::MatrixXd phiProjections(nlocal, nlocal);

          for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
          {
            Eigen::VectorXd phiProjectionsSingle = tripleProducts[basisIndex]*nIonVec;
            // Store components into rows of phiProjections matrix
            for (int colIndex = 0; colIndex < nlocal; colIndex++)
              phiProjections(basisIndex, colIndex) = phiProjectionsSingle(colIndex);
          }

          // construct arrays for passing into Petsc
          for (int k = 0; k < nlocal; k++)
            for (int m = 0; m < nlocal; m++)
              vals[nlocal*k + m] = phiProjections(k, m); // Default PetSc layout is row-major

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
          //for (unsigned k=0; k<nlocal; ++k)
          for (unsigned k=0; k<nlocal*nlocal; ++k)
            modVals[k] = 0.0;

// the following check is needed to ensure that the periodic BC
// // modifications to stiffness matrix are not applied twice
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

    // NOTE: This second loop is needed even though it is essentially the
    // same as the previous one as Petsc does not allow to call
    // MatZeroRows and MatSetValues without an intervening calls to
    // MatAssemblyBegin/MatAssemblyEnd. So if Petsc was not so annoying
    // this extra mess would not be needed. (Ammar, April 4 2012).
    for (int d = 0; d < NDIM; d++)
    {
      if (periodicFlgs[d] == true)
      {
        // fetch number of nodes on face of element
        unsigned nsl =  nodalBasis->getNumSurfUpperNodes(d);

        // space for mappings
        std::vector<int> lgSurfMap(nsl), lgLowerSurfMap(nsl);
        // space for local node numbers on faces
        std::vector<int> lgLocalNodeNum(nsl);

        // create region to loop over side
        Lucee::Region<NDIM, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
        // only update if we are on the correct ranks
        Lucee::Region<NDIM, int> defRgn = defRgnG; //.intersect(localRgn);

        // loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<NDIM> seqSide(defRgn);
        while (seqSide.step())
        {
          seqSide.fillWithIndex(idx);

          // set index into element basis
          nodalBasis->setIndex(idx);
          // get surface nodes -> global mapping
          nodalBasis->getSurfUpperLocalToGlobal(d, lgSurfMap);

          // reset corresponding rows (Note that some rows may be reset more
          // than once. This should not be a problem, though might make the
          // setup phase a bit slower).
#if PETSC_VERSION_GE(3,6,0)
          MatZeroRows(stiffMat, nsl, &lgSurfMap[0], 0.0, PETSC_NULL, PETSC_NULL);
#else
          MatZeroRows(stiffMat, nsl, &lgSurfMap[0], 0.0);
#endif

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

    // storage for computing source contribution
    std::vector<double> localMassSrc(nlocal);
    // create sequencer for looping over local box
    Lucee::RowMajorSequencer<NDIM> seqLocal(grid.getLocalRegion());

    // clear out existing stuff in source vector: this is required
    // otherwise successive calls to advance() will accumulate into source
    // from previous calls, which is of course not what we want.
    VecSet(globalSrc, 0.0);

    // loop, creating RHS (source terms)
    while (seqLocal.step())
    {
      seqLocal.fillWithIndex(idx);
      // set index into element basis
      nodalBasis->setIndex(idx);

      // get local mass matrix
      nodalBasis->getMassMatrix(localMassLucee);
      // copy to eigen
      copyLuceeToEigen(localMassLucee, massMatrix);

      // get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

      // Set inputs
      nElcIn.setPtr(nElcPtr, idx);
      nIonIn.setPtr(nIonPtr, idx);

      Eigen::VectorXd nDiffVec(nlocal);
      
      for (unsigned componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        nDiffVec(componentIndex) = nIonPtr[componentIndex] - nElcPtr[componentIndex];
        // Scale by various factors in the phi-equation
        nDiffVec(componentIndex) = Te0*nDiffVec(componentIndex)/(kPerpTimesRho*kPerpTimesRho);
      }

      // evaluate local mass matrix times local source
      Eigen::VectorXd rhsIntegrals = massMatrix*nDiffVec;

      // Copy over rhsIntegrals to localMassSrc
      for (int k = 0; k < nlocal; k++)
        localMassSrc[k] = rhsIntegrals(k);
      
      // accumulate it into Petsc vector
      VecSetValues(globalSrc, nlocal, &lgMap[0], &localMassSrc[0], ADD_VALUES);
    }

    // finish assembly of RHS
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

    // Now loop over each periodic direction, adding contribution from
    // last layer of cells on upper edges to lower edge nodes
    for (unsigned d = 0; d < NDIM; d++)
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
        Lucee::RowMajorSequencer<NDIM> seqSide(defRgn);
        while (seqSide.step())
        {
          seqSide.fillWithIndex(idx);

          // set index into element basis
          nodalBasis->setIndex(idx);
          // get local node number of upper face
          nodalBasis->getSurfUpperNodeNums(d, lgLocalNodeNum);

          // get local mass matrix
          nodalBasis->getMassMatrix(localMassLucee);
          // copy to eigen
          copyLuceeToEigen(localMassLucee, massMatrix);

          // Set inputs
          nElcIn.setPtr(nElcPtr, idx);
          nIonIn.setPtr(nIonPtr, idx);

          Eigen::VectorXd nDiffVec(nlocal);
          
          for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          {
            nDiffVec(componentIndex) = nIonPtr[componentIndex] - nElcPtr[componentIndex];
            // Scale by various factors in the phi-equation
            nDiffVec(componentIndex) = Te0*nDiffVec(componentIndex)/(kPerpTimesRho*kPerpTimesRho);
          }

          // evaluate local mass matrix times local source
          Eigen::VectorXd rhsIntegrals = massMatrix*nDiffVec;

          // Copy over rhsIntegrals to localMassSrc
          for (unsigned k=0; k<nlocal; ++k)
            localMassSrc[k] = rhsIntegrals(k);

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

    double resetVal[1];
    // reset source to apply boundary conditions
    std::map<int, double>::const_iterator rowItr
      = rowBcValues.begin();
    for ( ; rowItr != rowBcValues.end(); ++rowItr)
    {
      resetRow[0] = rowItr->first; // row index
      resetVal[0] = rowItr->second; // value
      VecSetValues(globalSrc, 1, resetRow, resetVal, INSERT_VALUES);
    }

    // reassemble RHS after application of Dirichlet Bcs
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

    // copy solution for use as initial guess in KSP solve
    copyFromGkeyllField(phiOut, initGuess);

#if PETSC_VERSION_GE(3,6,0)    
    KSPSetOperators(ksp, stiffMat, stiffMat);
#else
    KSPSetOperators(ksp, stiffMat, stiffMat, DIFFERENT_NONZERO_PATTERN);
#endif
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    
    KSPSetFromOptions(ksp);
    // now solve linear system (initGuess will contain solution)
    KSPSolve(ksp, globalSrc, initGuess);

    // check if solver converged
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    int itNum;
    KSPGetIterationNumber(ksp, &itNum);
    double resNorm;
    KSPGetResidualNorm(ksp, &resNorm);

    // construct message to send back to Lua
    std::ostringstream msgStrm;
    bool status = true;
    if (reason < 0) 
    {
      status = false;
      msgStrm << ElectrostaticContPhiUpdater::id << ": KSPSolve failed!";
      msgStrm << " Petsc reason code was " << reason << ".";
    }
    else
    {
      msgStrm << ElectrostaticContPhiUpdater::id << ": KSPSolve converged.";
    }
    msgStrm << " Number of iterations " << itNum
            << ". Final residual norm was " << resNorm;

    // copy solution from PetSc array to solution field
    copyFromPetscField(initGuess, phiOut);

    if (useCutoffVelocities == true)
    {
      // Dynvector containing the cutoff velocities computed at one or both edges
      const Lucee::DynVector<double>& cutoffVIn = this->getInp<Lucee::DynVector<double> >(2);
      // Compute cutoff velocity on the right edge
      std::vector<double> cutoffVelocities = cutoffVIn.getLastInsertedData();
      double phiS = 0.5*ELECTRON_MASS*cutoffVelocities[1]*cutoffVelocities[1]/ELEMENTARY_CHARGE;
    
      // Figure out how many exclusive nodes there are per cell
      std::vector<int> ndIds;
      nodalBasis->getExclusiveNodeIndices(ndIds);

      Lucee::Region<NDIM, int> globalRgnDup = phiOut.getExtRegion(); 
      // Add phiS to the solution
      for (int ix = globalRgnDup.getLower(0); ix < globalRgnDup.getUpper(0); ix++)
      {
        phiOut.setPtr(phiPtr, ix);

        for (int componentIndex = 0; componentIndex < ndIds.size(); componentIndex++)
          phiPtr[componentIndex] = phiPtr[componentIndex] + phiS;
      }
    }

    return Lucee::UpdaterStatus(status, std::numeric_limits<double>::max(),
      msgStrm.str());
  }

  void
  ElectrostaticContPhiUpdater::declareTypes()
  {
    // inputs n_e(x), n_i(x)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // optional input: cutoff velocity at edges
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // returns one output: phi(x) -- continuous version!
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  void
  ElectrostaticContPhiUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  void
  ElectrostaticContPhiUpdater::copyFromPetscField(Vec ptFld, Lucee::Field<NDIM, double>& gkFld)
  {
    PetscScalar *ptFldPtr;
#ifndef HAVE_MPI
      VecGetArray(ptFld, &ptFldPtr);
      nodalBasis->copyAllDataToField(ptFldPtr, gkFld);
      VecRestoreArray(ptFld, &ptFldPtr);
#else
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
#endif
  }

  void
  ElectrostaticContPhiUpdater::copyFromGkeyllField(const Lucee::Field<NDIM, double>& gkFld, Vec ptFld)
  {
    PetscScalar *ptFldPtr;
#ifndef HAVE_MPI
      VecGetArray(ptFld, &ptFldPtr);
      nodalBasis->copyAllDataFromField(gkFld, ptFldPtr);
      VecRestoreArray(ptFld, &ptFldPtr);
#else
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
// // layout is slightly different than Gkeyll layout. Most data is
// // already local and does not require communication.)
      VecScatterBegin(vecSctr, localData, ptFld, INSERT_VALUES, SCATTER_REVERSE);
      VecScatterEnd(vecSctr, localData, ptFld, INSERT_VALUES, SCATTER_REVERSE);
#endif
  }
}

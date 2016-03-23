/**
 * @file	LcElectromagneticContAUpdater.cpp
 *
 * @brief	Updater to compute A using continuous basis functions
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
#include <LcElectromagneticContAUpdater.h>
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

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *ElectromagneticContAUpdater::id = "ElectromagneticContAUpdater";

  ElectromagneticContAUpdater::ElectromagneticContAUpdater()
    : UpdaterIfc()
  {
  }

  ElectromagneticContAUpdater::~ElectromagneticContAUpdater()
  {
#ifdef PETSC_36    
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
  ElectromagneticContAUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerp"))
      kPerp = tbl.getNumber("kPerp");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify kPerp");

    if (tbl.hasNumber("elcMass"))
      elcMass = tbl.getNumber("elcMass");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify elcMass");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify ionMass");

    if (tbl.hasNumber("elcCharge"))
      elcCharge = tbl.getNumber("elcCharge");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify elcCharge");

    if (tbl.hasNumber("ionCharge"))
      ionCharge = tbl.getNumber("ionCharge");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify ionCharge");

    if (tbl.hasNumber("mu0"))
      mu0 = tbl.getNumber("mu0");
    else
      throw Lucee::Except("ElectromagneticContAUpdater::readInput: Must specify mu0");

    // force all directions to be periodic
    for (int i = 0; i < 1; i++) 
      periodicFlgs[i] = true;
  }

  void 
  ElectromagneticContAUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // global region to update
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<1> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[1];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();
    // number of global nodes
    int nglobal = nodalBasis->getNumGlobalNodes();

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
    throw Lucee::Except("ElectromagneticContAUpdater does not yet work in parallel!");
#else
    // Explicit initialization of stiffness matrix speeds up
    // initialization tremendously.
    int nz = 10; // number of non-zero entries per row (WHAT SHOULD IT REALLY BE?)
    if (nodalBasis->getNumNodes() > 4)
      nz = 21; // HACK FOR NOW
    MatCreateSeqAIJ(PETSC_COMM_SELF, nglobal, nglobal, nz, PETSC_NULL, &stiffMat);
#endif
    MatSetFromOptions(stiffMat);
    // Finalize assembly... is this needed?
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);
    
    // create and setup vector (NOTE: ACCORDING TO MIKE MCCOURT ONE SHOULD
    // USE MatGetVecs INSTEAD OF VecCreate. THIS ENSURES THAT THE PARALLEL
    // LAYOUT OF THE MATRIX & VECTOR ARE THE SAME)
    VecCreate(MPI_COMM_WORLD, &globalSrc);
    VecSetSizes(globalSrc, nglobal, PETSC_DECIDE);
    VecSetFromOptions(globalSrc);
    
    // finalize assembly
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);
    
    // create duplicate to store initial guess
    VecDuplicate(globalSrc, &initGuess); 

    // Create solver
    KSPCreate(MPI_COMM_WORLD, &ksp);
  }

  Lucee::UpdaterStatus 
  ElectromagneticContAUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& mom0ElcIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& mom0IonIn = this->getInp<Lucee::Field<1, double> >(1);
    const Lucee::Field<1, double>& mom1ElcIn = this->getInp<Lucee::Field<1, double> >(2);
    const Lucee::Field<1, double>& mom1IonIn = this->getInp<Lucee::Field<1, double> >(3);
    Lucee::Field<1, double>& aOut = this->getOut<Lucee::Field<1, double> >(0);

    Lucee::ConstFieldPtr<double> mom0ElcPtr = mom0ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom0IonPtr = mom0IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1ElcPtr = mom1ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1IonPtr = mom1IonIn.createConstPtr();
    Lucee::FieldPtr<double> aPtr = aOut.createPtr();

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
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion(); 
    // create sequencer for looping over local box
    Lucee::Region<1, int> localRgn = grid.getLocalRegion(); 
    Lucee::RowMajorSequencer<1> seq(localRgn);
    int idx[1];

    // loop, creating global stiffness matrix
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // set index into element basis
      nodalBasis->setIndex(idx);
      // Set inputs
      mom0ElcIn.setPtr(mom0ElcPtr, idx);
      mom0IonIn.setPtr(mom0IonPtr, idx);

      Eigen::VectorXd lhsVec(nlocal);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        lhsVec(componentIndex) = kPerp*kPerp + mu0*(
            elcCharge*elcCharge*mom0ElcPtr[componentIndex]/(elcMass*elcMass) + 
            ionCharge*ionCharge*mom0IonPtr[componentIndex]/(ionMass*ionMass) );

      Eigen::MatrixXd aProjections(nlocal, nlocal);

      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        Eigen::VectorXd aProjectionsSingle = tripleProducts[basisIndex]*lhsVec;
        // Store components into rows of aProjections matrix
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
          aProjections(basisIndex, colIndex) = aProjectionsSingle(colIndex);
      }

      // construct arrays for passing into Petsc
      for (int k = 0; k < nlocal; k++)
        for (int m = 0; m < nlocal; m++)
          vals[nlocal*k + m] = aProjections(k, m); // Default PetSc layout is row-major

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
    for (int d = 0; d < 1; d++)
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
        Lucee::Region<1, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
        // only update if we are on the correct ranks
        Lucee::Region<1, int> defRgn = defRgnG.intersect(localRgn);

        // loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<1> seqSide(defRgn);
        while (seqSide.step())
        {
          seqSide.fillWithIndex(idx);

          // set index into element basis
          nodalBasis->setIndex(idx);
          // Set inputs
          mom0ElcIn.setPtr(mom0ElcPtr, idx);
          mom0IonIn.setPtr(mom0IonPtr, idx);

          Eigen::VectorXd lhsVec(nlocal);

          for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
            lhsVec(componentIndex) = kPerp*kPerp + mu0*(
              elcCharge*elcCharge*mom0ElcPtr[componentIndex]/(elcMass*elcMass) + 
              ionCharge*ionCharge*mom0IonPtr[componentIndex]/(ionMass*ionMass) );

          Eigen::MatrixXd aProjections(nlocal, nlocal);

          for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
          {
            Eigen::VectorXd aProjectionsSingle = tripleProducts[basisIndex]*lhsVec;
            // Store components into rows of aProjections matrix
            for (int colIndex = 0; colIndex < nlocal; colIndex++)
              aProjections(basisIndex, colIndex) = aProjectionsSingle(colIndex);
          }

          // construct arrays for passing into Petsc
          for (int k = 0; k < nlocal; k++)
            for (int m = 0; m < nlocal; m++)
              vals[nlocal*k + m] = aProjections(k, m); // Default PetSc layout is row-major

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
          for (unsigned k=0; k<nlocal; ++k)
            modVals[k] = 0.0;

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
    for (int d = 0; d < 1; d++)
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
        Lucee::Region<1, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
        // only update if we are on the correct ranks
        Lucee::Region<1, int> defRgn = defRgnG.intersect(localRgn);

        // loop, modifying stiffness matrix
        Lucee::RowMajorSequencer<1> seqSide(defRgn);
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
#ifdef PETSC_36          
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

    // clear out existing stuff in source vector: this is required
    // otherwise successive calls to advance() will accumulate into source
    // from prevous calls, which is of course not what we want.
    VecSet(globalSrc, 0.0);

    // storage for computing source contribution
    std::vector<double> localMassSrc(nlocal);
    // create sequencer for looping over local box
    Lucee::RowMajorSequencer<1> seqLocal(grid.getLocalRegion());

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
      mom1ElcIn.setPtr(mom1ElcPtr, idx);
      mom1IonIn.setPtr(mom1IonPtr, idx);

      Eigen::VectorXd rhsVec(nlocal);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        rhsVec(componentIndex) = mu0*( 
            ionCharge*mom1IonPtr[componentIndex]/(ionMass*ionMass) + 
            elcCharge*mom1ElcPtr[componentIndex]/(elcMass*elcMass) );

      // evaluate local mass matrix times local source
      Eigen::VectorXd rhsIntegrals = massMatrix*rhsVec;

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
    for (int d = 0; d < 1; d++)
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
        Lucee::Region<1, int> defRgnG = 
          globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
        // only update if we are on the correct ranks
        Lucee::Region<1, int> defRgn = defRgnG.intersect(localRgn);

        // loop, modifying source vector
        Lucee::RowMajorSequencer<1> seqSide(defRgn);
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
          mom1ElcIn.setPtr(mom1ElcPtr, idx);
          mom1IonIn.setPtr(mom1IonPtr, idx);

          Eigen::VectorXd rhsVec(nlocal);

          for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
            rhsVec(componentIndex) = mu0*( 
                ionCharge*mom1IonPtr[componentIndex]/(ionMass*ionMass) + 
                elcCharge*mom1ElcPtr[componentIndex]/(elcMass*elcMass) );

          // evaluate local mass matrix times local source
          Eigen::VectorXd rhsIntegrals = massMatrix*rhsVec;

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
    PetscScalar *ptGuess;
    unsigned count = 0;
    VecGetArray(initGuess, &ptGuess);
    nodalBasis->copyAllDataFromField(aOut, ptGuess);
    VecRestoreArray(initGuess, &ptGuess);

#ifdef PETSC_36
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
      msgStrm << ElectromagneticContAUpdater::id << ": KSPSolve failed!";
      msgStrm << " Petsc reason code was " << reason << ".";
    }
    else
    {
      msgStrm << ElectromagneticContAUpdater::id << ": KSPSolve converged.";
    }
    msgStrm << " Number of iterations " << itNum
            << ". Final residual norm was " << resNorm;

    // copy solution from PetSc array to solution field
    PetscScalar *ptSol;
    VecGetArray(initGuess, &ptSol);
    nodalBasis->copyAllDataToField(ptSol, aOut);
    VecRestoreArray(initGuess, &ptSol);

    return Lucee::UpdaterStatus(status, std::numeric_limits<double>::max(),
      msgStrm.str());
  }

  void
  ElectromagneticContAUpdater::declareTypes()
  {
    // takes mom0Elc, mom0Ion, mom1Elc, mom1Ion (in p-space)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: A(x)
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ElectromagneticContAUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}

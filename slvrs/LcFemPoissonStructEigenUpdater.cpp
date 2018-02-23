/**
 * @file	LcFemPoissonStructEigenUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with FEM scheme on a structured grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFemPoissonStructEigenUpdater.h>
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
  using namespace Eigen;

  static const unsigned DX = 0;
  static const unsigned DY = 1;
  static const unsigned DZ = 2;
  static const unsigned LO = 0;
  static const unsigned HI = 1;
  static const unsigned NO_BC = 0;
  static const unsigned DIRICHLET_BC = 1;
  static const unsigned NEUMANN_BC = 2;
  static const unsigned DIRICHLET_VARIABLE_BC = 3;

  template <> const char *FemPoissonStructEigenUpdater<2>::id = "FemPoissonEigen2D";

  template <unsigned NDIM>
  FemPoissonStructEigenUpdater<NDIM>::FemPoissonStructEigenUpdater()
    : Lucee::UpdaterIfc()
  {
// this flag is needed to ensure arrays are only deleted if
// update() method is called at least once.
    runOnce = false;
  }


  template <unsigned NDIM>
  FemPoissonStructEigenUpdater<NDIM>::~FemPoissonStructEigenUpdater()
  {
    if (runOnce) 
    {
      //delete stiffMat;
      //delete globalSrc;
      //delete sourceModVec;
      //delete x;
      //delete solver;
    }
//    if (runOnce)
//    {
//#if PETSC_VERSION_GE(3,6,0)      
//      MatDestroy(&stiffMat);
//      VecDestroy(&globalSrc);
//      VecDestroy(&initGuess);
//      //VecDestroy(&localData);
//      //ISDestroy(&is);
//#else
//      MatDestroy(stiffMat);
//      VecDestroy(globalSrc);
//      VecDestroy(initGuess);
//      //VecDestroy(localData);
//      //ISDestroy(is);
//#endif      
//    }
  }

  template <unsigned NDIM>
  void
  FemPoissonStructEigenUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except(
        "FemPoissonStructEigenUpdater::readInput: Must specify element to use using 'basis'");

// check if source nodes are shared
    srcNodesShared = true;
    if (tbl.hasBool("sourceNodesShared"))
      srcNodesShared = tbl.getBool("sourceNodesShared");

// check if solution nodes are shared
    solNodesShared = true;
    if (tbl.hasBool("solutionNodesShared"))
      solNodesShared = tbl.getBool("solutionNodesShared");

// check if we are actually solving GK Poisson equations
    isGkPoisson = false;
    if (tbl.hasBool("isGyroKineticPoisson"))
      isGkPoisson = tbl.getBool("isGyroKineticPoisson");

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
          throw Lucee::Except("FemPoissonStructEigenUpdater:readInput: Incorrect 'periodicDirs'");
        }
        periodicFlgs[(unsigned) pd[i]] = true;
      }
    }

// check if all directions are periodic
    allPeriodic = true;
    for (unsigned i=0; i<NDIM; ++i)
      allPeriodic = allPeriodic && periodicFlgs[i];

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
    int nx = globalRgn.getShape(0);
    int ny = globalRgn.getShape(1);

// get BCs to apply
    if (tbl.hasTable("bcLeft"))
    {
      bc[DX][LO] = getBcData(tbl.getTable("bcLeft"));

      // get index of start and end of left boundary in global mapping
      // start
      nodalBasis->setIndex(0,1);
      bc[DX][LO].istart = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
      // end
      nodalBasis->setIndex(0,ny-1);
      bc[DX][LO].iend = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
    }
    if (tbl.hasTable("bcRight"))
    {
      bc[DX][HI] = getBcData(tbl.getTable("bcRight"));

      // get index of start and end of right boundary in global mapping
      // start
      nodalBasis->setIndex(nx,1);
      bc[DX][HI].istart = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
      // end
      nodalBasis->setIndex(nx,ny-1);
      bc[DX][HI].iend = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
    }
    if (NDIM>1)
    {
      if (tbl.hasTable("bcBottom"))
      {
        bc[DY][LO] = getBcData(tbl.getTable("bcBottom"));

        // get index of start and end of bottom boundary in global mapping
        // start
        nodalBasis->setIndex(0,0);
        bc[DY][LO].istart = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
        // end
        nodalBasis->setIndex(nx,0);
        bc[DY][LO].iend = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
      }
      if (tbl.hasTable("bcTop"))
      {
        bc[DY][HI] = getBcData(tbl.getTable("bcTop"));

        // get index of start and end of top boundary in global mapping
        // start
        nodalBasis->setIndex(0,ny);
        bc[DY][HI].istart = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
        // end
        nodalBasis->setIndex(nx,ny);
        bc[DY][HI].iend = nodalBasis->getLocalToGlobalInteriorBLRT(periodicFlgs);
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
          Lucee::Except lce("FemPoissonStructEigenUpdater:: Must specify BCs on each side");
          throw lce;
        }
      }
      else
      {
// ensure no BCs are specified
        if (bc[dir][0].isSet || bc[dir][1].isSet)
        {
          Lucee::Except lce("FemPoissonStructEigenUpdater:: Cannot specify BCs");
          lce << " if a direction is periodic";
          throw lce;
        }
      }
    }

    modifierConstant = 0.0;
// read in modifier constant (if any)
    if (tbl.hasNumber("modifierConstant"))
      modifierConstant = tbl.getNumber("modifierConstant");

    laplacianWeight = 1.0;
// read in laplacian weight (if any)
    if (tbl.hasNumber("laplacianWeight"))
      laplacianWeight = tbl.getNumber("laplacianWeight");

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
  FemPoissonStructEigenUpdater<NDIM>::initialize()
  {
//#define DMSG(s) std::cout << "Rank " << this->getComm()->getRank() << ": " << s << std::endl;
#define DMSG(s) ;

    DMSG("Inside initialize");

    Lucee::UpdaterIfc::initialize();

    int nglobal = nodalBasis->getNumGlobalNodes(periodicFlgs);
    unsigned nlocal = nodalBasis->getNumNodes();
    int ninterior = nodalBasis->getNumInteriorNodes();
    int nz = nlocal*(std::pow(2.0, 1.0*NDIM)+1); // number of nonzeros per row

    Lucee::Matrix<double> localStiff(nlocal, nlocal), localMass(nlocal, nlocal);
    std::vector<int> lgMap(nlocal);
    std::vector<Triplet<double> > tripletList, identityTripletList;
    tripletList.reserve(nz*nglobal); // estimate number of nonzero elements
    DMSG("Matrix and vectors allocated");

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion(); 
// sequence over global region, so that each proc builds entire global stiffness matrix
    Lucee::RowMajorSequencer<NDIM> seq(globalRgn); 
    int idx[NDIM];

// loop, creating global stiffness matrix
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);

      if (isGkPoisson)
        nodalBasis->getPerpStiffnessMatrix(localStiff);
      else
        nodalBasis->getStiffnessMatrix(localStiff);
      if(modifierConstant!=0.) nodalBasis->getMassMatrix(localMass);

// get local-to-global interior-boundary mapping
// this takes care of periodic boundary conditions by mapping 
// e.g. right boundary to left boundary for periodicity in x
      nodalBasis->getLocalToGlobalInteriorBLRT(lgMap, periodicFlgs);
// make triplets for constructing Eigen SparseMatrix
      for (unsigned k=0; k<nlocal; ++k)
      {
        for (unsigned m=0; m<nlocal; ++m) {
          double val = -laplacianWeight*localStiff(k,m);
          if(modifierConstant!=0.) val += modifierConstant*localMass(k,m);
          unsigned globalIdx_k = lgMap[k];
          unsigned globalIdx_m = lgMap[m];
          
          tripletList.push_back(Triplet<double>(globalIdx_k, globalIdx_m, val)); 
        }
      }
    }

// construct SparseMatrix from triplets
    SparseMatrix<double,RowMajor> stiffMatRowMajor(nglobal, nglobal);
    stiffMatRowMajor.setFromTriplets(tripletList.begin(), tripletList.end());

    DMSG("Beginning modification for Dirichlet BCs");
// handle Dirichlet BCs
// note: we do not need to do anything for Neumann BCs

    int nDirichlet = 0;
// set rows corresponding to Dirichlet BCs to zero
    for (unsigned d=0; d<NDIM; ++d)
    {
      for (unsigned side=0; side<2; ++side)
      {
        if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC ||
            bc[d][side].isSet && bc[d][side].type == DIRICHLET_VARIABLE_BC)
        {
          unsigned start = bc[d][side].istart;
          unsigned end = bc[d][side].iend;
          int nRows = end-start;
          stiffMatRowMajor.middleRows(start, nRows) = SparseMatrix<double,RowMajor>(nRows, stiffMatRowMajor.cols());   

          nDirichlet += nRows;
        }
      }
    }
// if all periodic, make the bottom left corner effectively Dirichlet
    if(adjustSource) 
    {
      stiffMatRowMajor.middleRows(ninterior, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
    }

// create column major copy of stiffMat so that we can zero columns
    stiffMat = SparseMatrix<double,ColMajor>(stiffMatRowMajor);
// de-allocate row major copy
    stiffMatRowMajor.resize(0,0);

// create matrix for modifying source when there are Dirichlet BCs
// this matrix will have nonzero entries only for the blocks consisting
// of Dirichlet cols and non-Dirichlet rows
    SparseMatrix<double,ColMajor> sourceModMat(nglobal, nglobal);

// sparse matrix for identity blocks in rows/cols corresponding to Dirichlet BCs
    SparseMatrix<double,ColMajor> dirichletIdentity(nglobal, nglobal);
    identityTripletList.reserve(nDirichlet); // estimate number of nonzero elements

    for (unsigned d=0; d<NDIM; ++d)
    {
      for (unsigned side=0; side<2; ++side)
      {
        if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC ||
            bc[d][side].isSet && bc[d][side].type == DIRICHLET_VARIABLE_BC)
        {
          unsigned start = bc[d][side].istart;
          unsigned end = bc[d][side].iend;
          int nCols = end-start;
// copy cols corresponding to Dirichlet BCs to sourceModMat
          sourceModMat.middleCols(start, nCols) = stiffMat.middleCols(start, nCols);

// in stiffMat, set cols corresponding to Dirichlet BCs to zero
          stiffMat.middleCols(start, nCols) = SparseMatrix<double,ColMajor>(stiffMat.rows(), nCols);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
          for(unsigned i=start; i<start+nCols; i++) {
            identityTripletList.push_back(Triplet<double>(i, i, 1.)); 
          } 
        }
      }
    }
// if all periodic, make the bottom left corner effectively Dirichlet
    if(adjustSource) 
    {
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(ninterior, 1) = stiffMat.middleCols(ninterior, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(ninterior, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(ninterior, ninterior, 1.));
    }
    dirichletIdentity.setFromTriplets(identityTripletList.begin(), identityTripletList.end());
    stiffMat+=dirichletIdentity;
    
// get rid of zeros in sparse matrices
    stiffMat.prune(0.);
    sourceModMat.prune(0.);

// create vector of Dirichlet values
    SparseVector<double> dirichletVec(nglobal);
    dirichletVec.reserve(nDirichlet+1);
    for (unsigned d=0; d<NDIM; ++d)
    {
      for (unsigned side=0; side<2; ++side)
      {
        if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC)
        {
          for(unsigned i=bc[d][side].istart; i<bc[d][side].iend; i++) {
            dirichletVec.coeffRef(i) = bc[d][side].value;
          }
        }
      }
    }
    if(adjustSource)
    {
      // if all periodic, set bottom left corner to 0
      dirichletVec.coeffRef(ninterior) = 0.;
    }

// calculate vector to subtract from source
    sourceModVec = SparseVector<double>(nglobal);
    sourceModVec = sourceModMat*dirichletVec;

    DMSG("Modification for Dirichlet BCs completed");

    
    if (writeMatrix)
    {
      std::string outName = Loki::SingletonHolder<Lucee::Globals>
        ::Instance().outPrefix + "-poisson-stiffnessMatrix";
      saveMarket(stiffMat, outName, Symmetric);
    }

// initialize solver
    //solver = SimplicialLDLT<SparseMatrix<double> >;

// compute step: reorder and factorize stiffMat to prepare for solve 
    solver.compute(stiffMat);

// initialize global source vector
    globalSrc = VectorXd(nglobal);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FemPoissonStructEigenUpdater<NDIM>::update(double t)
  {
// set flag to indicate we should delete PetSc vectors/matrices
    runOnce = true;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& src = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& sol = this->getOut<Lucee::Field<NDIM, double> >(0);

 // NEED TO MAKE SURE ENTIRE SRC IS ON ALL PROCS

    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion(); 
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned nglobal = nodalBasis->getNumGlobalNodes(periodicFlgs);
    int ninterior = nodalBasis->getNumInteriorNodes();

    Lucee::Matrix<double> localMass(nlocal, nlocal);
    std::vector<int> lgMap(nlocal);
    std::vector<double> localSrc(nlocal), localMassSrc(nlocal);
    Lucee::ConstFieldPtr<double> srcPtr = src.createConstPtr();
    Lucee::FieldPtr<double> solPtr = sol.createPtr();

// every proc does its own global solve
    Lucee::RowMajorSequencer<NDIM> seq(grid.getGlobalRegion());
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
    globalSrc.setZero(nglobal);

// loop, creating RHS (source terms)
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);

      nodalBasis->getMassMatrix(localMass);
      nodalBasis->getLocalToGlobalInteriorBLRT(lgMap, periodicFlgs);

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
// map into global source vector
      for(unsigned i=0; i<nlocal; i++)
      {
        globalSrc.coeffRef(lgMap[i]) += localMassSrc[i];
      }
    }

// replace dirichlet nodes of global source with dirichlet values
    for (unsigned d=0; d<NDIM; ++d)
    {
      for (unsigned side=0; side<2; ++side)
      {
        if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC)
        {
          for(unsigned i=bc[d][side].istart; i<bc[d][side].iend; i++) {
            globalSrc.coeffRef(i) = bc[d][side].value;
          }
        }
      }
    }
    if(adjustSource) {
      globalSrc.coeffRef(ninterior) = 0.;
    }
 
// modify non-dirichlet rows of source by subtracting sourceModVec
    globalSrc -= sourceModVec;

    if (writeMatrix)
    {
      std::string outName = Loki::SingletonHolder<Lucee::Globals>
        ::Instance().outPrefix + "-poisson-globalSrc";
      saveMarketVector(globalSrc, outName);
    }

// solve linear system
    x = solver.solve(globalSrc);

    if (writeMatrix)
    {
      std::string outName = Loki::SingletonHolder<Lucee::Globals>
        ::Instance().outPrefix + "-poisson-solution";
      saveMarketVector(x, outName);
    }

// remap global solution to local gkeyll solution field
// only need to sequence over local proc region
    Lucee::RowMajorSequencer<NDIM> localseq(grid.getLocalRegion());
    while (localseq.step())
    {
      localseq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);

      nodalBasis->getLocalToGlobalInteriorBLRT(lgMap, periodicFlgs);

      sol.setPtr(solPtr, idx);
      for (unsigned k=0; k<nlocal; ++k) {
        solPtr[k] = x.coeffRef(lgMap[k]);
      }
    }

    bool status = true;
    return Lucee::UpdaterStatus(status, std::numeric_limits<double>::max(),
      "");
  }

  template <unsigned NDIM>
  void
  FemPoissonStructEigenUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>)); // optional, for variable dirichlet
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  typename FemPoissonStructEigenUpdater<NDIM>::FemPoissonBcData
  FemPoissonStructEigenUpdater<NDIM>::getBcData(const Lucee::LuaTable& bct) const
  {
    FemPoissonBcData bcData;
    if (bct.getString("T") == "D")
      bcData.type = DIRICHLET_BC;
    else if ((bct.getString("T") == "N"))
      bcData.type = NEUMANN_BC;
    else if (bct.getString("T") == "D_VAR")
      bcData.type = DIRICHLET_VARIABLE_BC;
    else
    {
      Lucee::Except lce(
        "FemPoissonStructEigenUpdater::readInput: Must specify one of \"D\", \"D_VAR\", or \"N\". ");
      lce << "Specified \"" << bct.getString("T") << " instead";
      throw lce;
    }
    bcData.value = bct.getNumber("V");
    bcData.isSet = true;
    
    return bcData;
  }

  template <unsigned NDIM>
  double
  FemPoissonStructEigenUpdater<NDIM>::getFieldIntegral(const Lucee::Field<NDIM, double>& fld, 
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

/*
  template <unsigned NDIM>
  void
  FemPoissonStructEigenUpdater<NDIM>::copyFromPetscField(Vec ptFld, Lucee::Field<NDIM, double>& gkFld)
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
  FemPoissonStructEigenUpdater<NDIM>::copyFromGkeyllField(const Lucee::Field<NDIM, double>& gkFld, Vec ptFld)
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
*/

// instantiations
  template class FemPoissonStructEigenUpdater<1>;
  template class FemPoissonStructEigenUpdater<2>;
  template class FemPoissonStructEigenUpdater<3>;
}

// The only place solution appears is when setting an initial
// condition for iterative solver and when copying solution back to
// Gkeyll fields. This is now the remaining problem in getting the
// correct solution in parallel.
//
// One thing I don't understand is how the local->global mapping can
// work.

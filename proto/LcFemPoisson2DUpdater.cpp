/**
 * @file	LcFemPoisson2DUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with FEM scheme in 1D.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFemPoisson2DUpdater.h>
#include <LcMatrix.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *FemPoisson2DUpdater::id = "FemPoisson2D";

  FemPoisson2DUpdater::FemPoisson2DUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  FemPoisson2DUpdater::~FemPoisson2DUpdater()
  { // get rid of stuff
    MatDestroy(stiffMat);
    VecDestroy(globalSrc);
    VecDestroy(initGuess);
  }

  void
  FemPoisson2DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc>("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc>("basis");
    else
      throw Lucee::Except("FemPoisson2DUpdater::readInput: Must specify element to use using 'basis'");
  }

  void
  FemPoisson2DUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// number of global nodes
    unsigned nglobal = nodalBasis->getNumGlobalNodes();
// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// now create and initialize Petsc matrix to store stiffness matrix (LHS for Poisson equation)
    MatCreate(MPI_COMM_WORLD, &stiffMat);
    MatSetSizes(stiffMat, nglobal, nglobal, PETSC_DECIDE, PETSC_DECIDE);
    MatSetFromOptions(stiffMat);

// create and setup vector
    VecCreate(MPI_COMM_WORLD, &globalSrc);
    VecSetSizes(globalSrc, nglobal, PETSC_DECIDE);
    VecSetFromOptions(globalSrc);

// get global region
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

// local stiffness matrix
    Lucee::Matrix<double> localStiff(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// storage for passing to petsc
    std::vector<PetscScalar> vals(nlocal*nlocal);

// loop, creating stiffness matrix
    for (int i=globalRgn.getLower(0); i<globalRgn.getUpper(0); ++i)
    {
      nodalBasis->setIndex(i);

// get local stiffness matrix
      nodalBasis->getStiffnessMatrix(localStiff);
// get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

// construct arrays for passing into Petsc
      for (unsigned k=0; k<nlocal; ++k)
      {
        for (unsigned m=0; m<nlocal; ++m)
          vals[nlocal*k+m] = -localStiff(k,m);
      }

// insert into global stiffness matrix, adding them to existing value
      MatSetValues(stiffMat, nlocal, &lgMap[0], nlocal, &lgMap[0], 
        &vals[0], ADD_VALUES);
    }

    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

//  finalize vector and matrix assembly
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

// create KSP context
    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, stiffMat, stiffMat, DIFFERENT_NONZERO_PATTERN);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetFromOptions(ksp);

// UNCOMMENT FOLLOWING LINE TO VIEW STIFFNESS MATRIX
//    MatView(stiffMat, PETSC_VIEWER_STDOUT_SELF);

// create duplicate to store initial guess
    VecDuplicate(globalSrc, &initGuess);
  }

  Lucee::UpdaterStatus
  FemPoisson2DUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input/output fields
    const Lucee::Field<2, double>& src = this->getInp<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& sol = this->getOut<Lucee::Field<2, double> >(0);

// number of global nodes
    unsigned nglobal = nodalBasis->getNumGlobalNodes();
// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// global region to index
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

// local stiffness matrix
    Lucee::Matrix<double> localMass(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// storage for computing source contribution
    std::vector<double> localSrc(nlocal), localMassSrc(nlocal);

    Lucee::ConstFieldPtr<double> srcPtr = src.createConstPtr();
    Lucee::ConstFieldPtr<double> srcPtrp = src.createConstPtr();
// loop, creating RHS (source terms)
    for (int i=globalRgn.getLower(0); i<globalRgn.getUpper(0); ++i)
    {
      nodalBasis->setIndex(i);
      src.setPtr(srcPtr, i);
      src.setPtr(srcPtrp, i+1);

// get local mass matrix
      nodalBasis->getMassMatrix(localMass);
// get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

// now compute source at each local node (cell i does not own the right-most node)
      for (unsigned k=0; k<nlocal-1; ++k)
        localSrc[k] = srcPtr[k];
      localSrc[nlocal-1] = srcPtrp[0]; // get right-most data from first node of right cell

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

    int firstLast[2];
    firstLast[0] = 0; firstLast[1] = nglobal-1;
    double vals[2];
    vals[0] = leftEdge; vals[1] = rightEdge;
// replace first/last values with specified BC values
    VecSetValues(globalSrc, 2, firstLast, vals, INSERT_VALUES);

// NOTE: We need to reassemble otherwise PetSc barfs. THIS PROBABLY IS
// NOT THE BEST WAY TO DO THIS STUFF, ANYWAY (Ammar, 3/07/2012)
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

// UNCOMMENT FOLLOWING LINE TO VIEW RHS VECTOR
//    VecView(globalSrc, PETSC_VIEWER_STDOUT_SELF);

// copy solution for use as initial guess in KSP solve
    Lucee::ConstFieldPtr<double> solConstPtr = sol.createConstPtr();
    PetscScalar *ptGuess;
    unsigned count = 0;
    VecGetArray(initGuess, &ptGuess);

    for (int i=globalRgn.getLower(0); i<globalRgn.getUpper(0); ++i)
    {
      sol.setPtr(solConstPtr, i);
      for (unsigned k=0; k<nlocal-1; ++k)
        ptGuess[count++] = solConstPtr[k];
    }
// copy solution at last node
    sol.setPtr(solConstPtr, globalRgn.getUpper(0));
    ptGuess[count] = solConstPtr[0];

    VecRestoreArray(initGuess, &ptGuess);

// now solve linear system (initGuess will contain solution)
    KSPSolve(ksp, globalSrc, initGuess);

    Lucee::FieldPtr<double> solPtr = sol.createPtr();
    count = 0;
    PetscScalar *ptSol;
    VecGetArray(initGuess, &ptSol);

    for (int i=globalRgn.getLower(0); i<globalRgn.getUpper(0); ++i)
    {
      sol.setPtr(solPtr, i);
      for (unsigned k=0; k<nlocal-1; ++k)
        solPtr[k] = ptSol[count++];
    }
// copy solution at last node
    sol.setPtr(solPtr, globalRgn.getUpper(0));
    solPtr[0] = ptSol[count];

    VecRestoreArray(initGuess, &ptSol);

    return Lucee::UpdaterStatus();
  }

  void
  FemPoisson2DUpdater::declareTypes()
  {
// takes one input (source terms)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// returns one output, solution
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}

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
  static const unsigned DX = 0;
  static const unsigned DY = 1;
  static const unsigned DZ = 2;
  static const unsigned LO = 0;
  static const unsigned HI = 1;
  static const unsigned DIRICHLET_BC = 0;
  static const unsigned NEUMANN_BC = 1;

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
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("FemPoisson2DUpdater::readInput: Must specify element to use using 'basis'");

// get BCs to apply
    if (tbl.hasTable("bcLeft"))
    {
      bc[DX][LO] = getBcData(tbl.getTable("bcLeft"));
    }
    if (tbl.hasTable("bcRight"))
    {
      bc[DX][HI] = getBcData(tbl.getTable("bcRight"));
    }
    if (tbl.hasTable("bcBottom"))
    {
      bc[DY][LO] = getBcData(tbl.getTable("bcBottom"));
    }
    if (tbl.hasTable("bcTop"))
    {
      bc[DY][HI] = getBcData(tbl.getTable("bcTop"));
    }

// TODO: NEED TO CORRECT THESE TESTS
// some sanity checks: must either specify both BCs along a direction
// or none. In the latter case the BC is assumed to be periodic.

//     if (bcX.size() > 0 && bcX.size() != 2)
//       throw Lucee::Except(
//         "FemPoisson2DUpdater::readInput: Must specify both bcLeft/bcRight");
//     if (bcY.size() > 0 && bcY.size() != 2)
//       throw Lucee::Except(
//         "FemPoisson2DUpdater::readInput: Must specify both bcBottom/bcTop");
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

// local stiffness matrix
    Lucee::Matrix<double> localStiff(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// storage for passing to petsc
    std::vector<PetscScalar> vals(nlocal*nlocal);

    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// create sequencer for looping over local box
    Lucee::Region<2, int> localRgn = grid.getLocalRegion(); 
    Lucee::RowMajorSequencer<2> seq(localRgn);
    int idx[2];

// loop, creating stiffness matrix
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// set index into element basis
      nodalBasis->setIndex(idx);

// get local stiffness matrix
      nodalBasis->getStiffnessMatrix(localStiff);
// get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

// construct arrays for passing into Petsc
      for (unsigned k=0; k<nlocal; ++k)
      {
        for (unsigned m=0; m<nlocal; ++m)
          vals[nlocal*k+m] = -localStiff(k,m); // Default PetSc layout is row-major
      }

// insert into global stiffness matrix, adding them to existing value
      MatSetValues(stiffMat, nlocal, &lgMap[0], nlocal, &lgMap[0],
        &vals[0], ADD_VALUES);
    }

    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

// now modify values in stiffness matrix based on Dirichlet Bcs
    for (unsigned d=0; d<2; ++d)
    { // d ranges over dimension
      for (unsigned side=0; side<2; ++side)
      { // side ranges over lower/upper side along d
        if (bc[d][side].type == NEUMANN_BC)
          break; // do nothing for Neumann Bcs

// fetch number of nodes on face of element
        unsigned nsl = side==0 ? 
          nodalBasis->getNumSurfLowerNodes(d) : nodalBasis->getNumSurfUpperNodes(d);

// allocate space for mapping
        std::vector<int> lgSurfMap(nsl);

        double dv = bc[d][side].value;
// create region to loop over side
        Lucee::Region<2, int> defRgn = localRgn.deflate(d);
        Lucee::RowMajorSequencer<2> seq(defRgn);
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
        }
      }
    }

// reassemble matrix after modification
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

// local mass matrix
    Lucee::Matrix<double> localMass(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// storage for computing source contribution
    std::vector<double> localSrc(nlocal), localMassSrc(nlocal);

// create sequencer for looping over local box
    Lucee::RowMajorSequencer<2> seq(grid.getLocalRegion());
    int idx[2];
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

// now compute source at each local node
      nodalBasis->extractFromField(src, localSrc);

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

// UNCOMMENT FOLLOWING LINE TO VIEW RHS VECTOR
//    VecView(globalSrc, PETSC_VIEWER_STDOUT_SELF);

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
    bool status = true;
    if (reason < 0) status = false;

// copy solution from PetSc array to solution field
    PetscScalar *ptSol;
    VecGetArray(initGuess, &ptSol);
    nodalBasis->copyAllDataToField(ptSol, sol);
    VecRestoreArray(initGuess, &ptSol);

    return Lucee::UpdaterStatus(status, std::numeric_limits<double>::max());
  }

  void
  FemPoisson2DUpdater::declareTypes()
  {
// takes one input (source terms)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// returns one output, solution
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  FemPoisson2DUpdater::FemPoissonBcData
  FemPoisson2DUpdater::getBcData(const Lucee::LuaTable& bct) const
  {
    FemPoissonBcData bcData;
    if (bct.getString("T") == "D")
      bcData.type = DIRICHLET_BC;
    else if ((bct.getString("T") == "N"))
      bcData.type = NEUMANN_BC;
    else
    {
      Lucee::Except lce(
        "FemPoisson2DUpdater::readInput: Must specify one of \"D\" or \"N\". ");
      lce << "Specified \"" << bct.getString("T") << " instead";
      throw lce;
    }
    return bcData;
  }
}

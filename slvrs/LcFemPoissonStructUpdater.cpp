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
      throw Lucee::Except("FemPoissonStructUpdater::readInput: Must specify element to use using 'basis'");

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

// TODO: NEED TO CORRECT THESE TESTS
// some sanity checks: must either specify both BCs along a direction
// or none. In the latter case the BC is assumed to be periodic.

//     if (bcX.size() > 0 && bcX.size() != 2)
//       throw Lucee::Except(
//         "FemPoissonStructUpdater::readInput: Must specify both bcLeft/bcRight");
//     if (bcY.size() > 0 && bcY.size() != 2)
//       throw Lucee::Except(
//         "FemPoissonStructUpdater::readInput: Must specify both bcBottom/bcTop");
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

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
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
    for (unsigned d=0; d<NDIM; ++d)
    {
      for (unsigned side=0; side<2; ++side)
      {
        if (bc[d][side].type == NEUMANN_BC)
          break; // do nothing for Neumann Bcs

// fetch number of nodes on face of element
        unsigned nsl = side==0 ? 
          nodalBasis->getNumSurfLowerNodes(d) : nodalBasis->getNumSurfUpperNodes(d);

// allocate space for mapping
        std::vector<int> lgSurfMap(nsl);

        double dv = bc[d][side].value;
// create region to loop over side
        Lucee::Region<NDIM, int> defRgn = side==0 ? 
          localRgn.resetBounds(d, localRgn.getLower(d), localRgn.getLower(d)+1) :
          localRgn.resetBounds(d, localRgn.getUpper(d)-1, localRgn.getUpper(d)) ;

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

// reassemble matrix after modification
    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);

// UNCOMMENT FOLLOWING LINE TO VIEW STIFFNESS MATRIX
//    MatView(stiffMat, PETSC_VIEWER_STDOUT_SELF);

//  finalize vector and matrix assembly
    VecAssemblyBegin(globalSrc);
    VecAssemblyEnd(globalSrc);

// create KSP context
    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, stiffMat, stiffMat, DIFFERENT_NONZERO_PATTERN);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetFromOptions(ksp);

// create duplicate to store initial guess
    VecDuplicate(globalSrc, &initGuess);
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
    Lucee::RowMajorSequencer<NDIM> seq(grid.getLocalRegion());
    int idx[NDIM];
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

    int resetRow[1];
    double resetVal[1];
// reset source to apply Dirichlet boundary conditions
    std::map<int, double>::const_iterator rItr
      = rowBcValues.begin();
    for ( ; rItr != rowBcValues.end(); ++rItr)
    {
      resetRow[0] = rItr->first; // row index
      resetVal[0] = rItr->second; // Dirichlet value
      VecSetValues(globalSrc, 1, resetRow, resetVal, INSERT_VALUES);
    }

// reassemble RHS after application of Dirichlet Bcs
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
    
    return bcData;
  }

// instantiations
  template class FemPoissonStructUpdater<1>;
  template class FemPoissonStructUpdater<2>;
  template class FemPoissonStructUpdater<3>;
}

/**
 * @file	LcFemPoisson1DUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with FEM scheme in 1D.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFemPoisson1DUpdater.h>
#include <LcMatrix.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *FemPoisson1DUpdater::id = "FemPoisson1D";

  FemPoisson1DUpdater::FemPoisson1DUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  void
  FemPoisson1DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc>("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc>("basis");
    else
      throw Lucee::Except("FemPoisson1DUpdater::readInput: Must specify element to use using 'basis'");

// get potentials on left and right edges
    leftEdge = tbl.getNumber("leftEdge");
    rightEdge = tbl.getNumber("rightEdge");
  }

  void
  FemPoisson1DUpdater::initialize()
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
    VecCreate(MPI_COMM_WORLD, &phin);
    VecSetSizes(phin, nglobal, PETSC_DECIDE);
    VecSetFromOptions(phin);

// get global region
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

// local stiffness matrix
    Lucee::Matrix<double> localStiff(nlocal, nlocal);
// map for local indices to global indices
    std::vector<int> lgMap(nlocal);

// loop, creating stiffness matrix
    for (int i=globalRgn.getLower(0); i<globalRgn.getUpper(1); ++i)
    {
      nodalBasis->setIndex(i);

// get local stiffness matrix
      nodalBasis->getStiffnessMatrix(localStiff);
// get local to global mapping
      nodalBasis->getLocalToGlobal(lgMap);

// insert values into global stiffness matrix
    }

//  finalize vector and matrix assembly
    VecAssemblyBegin(phin);
    VecAssemblyEnd(phin);

    MatAssemblyBegin(stiffMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiffMat, MAT_FINAL_ASSEMBLY);
    
// create duplicates
    VecDuplicate(phin, &rhs);
    MatDuplicate(stiffMat, MAT_COPY_VALUES, &lhs);
  }

  Lucee::UpdaterStatus
  FemPoisson1DUpdater::update(double t)
  {

    return Lucee::UpdaterStatus();
  }

  void
  FemPoisson1DUpdater::declareTypes()
  {
// takes one input (source terms)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
// returns one output, solution
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}

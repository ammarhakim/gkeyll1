/**
 * @file	LcNodalDisContHyperUpdater.cpp
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcDirSequencer.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcMathLib.h>
#include <LcNodalDisContHyperUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *NodalDisContHyperUpdater<1>::id = "NodalDgHyper1D";
  template <> const char *NodalDisContHyperUpdater<2>::id = "NodalDgHyper2D";
  template <> const char *NodalDisContHyperUpdater<3>::id = "NodalDgHyper3D";

  template <unsigned NDIM>
  NodalDisContHyperUpdater<NDIM>::NodalDisContHyperUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>  
  void 
  NodalDisContHyperUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDisContHyperUpdater::readInput: Must specify element to use using 'basis'");

// equation to solve
    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("NodalDisContHyperUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around
  }

  template <unsigned NDIM>
  void 
  NodalDisContHyperUpdater<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx); // fetch index of first location
    nodalBasis->setIndex(idx); // set index into basis function
    
    unsigned nlocal = nodalBasis->getNumNodes();

// get node numbers on each lower and upper edges
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      lowerNodeNums[dir].nums.resize(nodalBasis->getNumSurfLowerNodes(dir));
      nodalBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

// reset numbers as element offsets them with 1
      for (unsigned k=0; k<nodalBasis->getNumSurfLowerNodes(dir); ++k)
        lowerNodeNums[dir].nums[k] += -1;

      upperNodeNums[dir].nums.resize(nodalBasis->getNumSurfUpperNodes(dir));
      nodalBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);

// reset numbers as element offsets them with 1
      for (unsigned k=0; k<nodalBasis->getNumSurfUpperNodes(dir); ++k)
        upperNodeNums[dir].nums[k] += -1;
    }

// space for mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// get stiffness matrix
      stiffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, stiffMatrix[dir].m);

// multiply matrices by inverse of mass matrix
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, stiffMatrix[dir].m);

// compute lift matrices
      lowerLift[dir].m = Lucee::Matrix<double>(nlocal, 
        nodalBasis->getNumSurfLowerNodes(dir));

      nodalBasis->getLowerFaceMassMatrix(dir, lowerLift[dir].m);
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, lowerLift[dir].m);
      
      upperLift[dir].m = Lucee::Matrix<double>(nlocal, 
        nodalBasis->getNumSurfUpperNodes(dir));

      nodalBasis->getUpperFaceMassMatrix(dir, upperLift[dir].m);
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, upperLift[dir].m);
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  NodalDisContHyperUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// get input field
    const Lucee::Field<NDIM, double>& q = this->getInp<Lucee::Field<NDIM, double> >(0);
// get output field
    Lucee::Field<NDIM, double>& qNew = this->getOut<Lucee::Field<NDIM, double> >(0);

// time-step
    double dt = t-this->getCurrTime();
// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

// maximum CFL number used
    double cfla = 0.0;

    unsigned nlocal = nodalBasis->getNumNodes(); // local nodes
    unsigned meqn = equation->getNumEqns();// num of equations

    qNew = 0.0; // clear out contents

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
// compute contribution from volume integrals
    while (seq.step())
    {
    }

// compute contributions from surface integrals

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>  
  void
  NodalDisContHyperUpdater<NDIM>::declareTypes()
  {
// takes one input
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// returns one output
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void 
  NodalDisContHyperUpdater<NDIM>::matVec(double m, const Lucee::Matrix<double>& mat,
    const double* vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned i=0; i<rows; ++i)
    {
      tv = 0.0;
      for (unsigned j=0; j<cols; ++j)
        tv += mat(i,j)*vec[j];
      out[i] = m*tv + v*out[i];
    }
  }

// instantiations
  template class NodalDisContHyperUpdater<1>;
  template class NodalDisContHyperUpdater<2>;
  template class NodalDisContHyperUpdater<3>;
}

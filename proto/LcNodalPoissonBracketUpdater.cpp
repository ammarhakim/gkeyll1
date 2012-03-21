/**
 * @file	LcNodalPoissonBracketUpdater.cpp
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinAlgebra.h>
#include <LcNodalPoissonBracketUpdater.h>

namespace Lucee
{
  const char *NodalPoissonBracketUpdater::id = "PoissonBracket";

  NodalPoissonBracketUpdater::NodalPoissonBracketUpdater()
    : UpdaterIfc(), diffMatrix_x(1,1), diffMatrix_y(1,1),
      stiffMatrix_x(1,1), stiffMatrix_y(1,1)
  {
  }
  
  void 
  NodalPoissonBracketUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify element to use using 'basis'");
  }

  void 
  NodalPoissonBracketUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(localRgn.getLower(0), localRgn.getLower(1));

// compute various matrices
    unsigned nlocal = nodalBasis->getNumNodes();

    stiffMatrix_x = Lucee::Matrix<double>(nlocal, nlocal);
    stiffMatrix_y = Lucee::Matrix<double>(nlocal, nlocal);

    nodalBasis->getGradStiffnessMatrix(0, stiffMatrix_x);
    nodalBasis->getGradStiffnessMatrix(1, stiffMatrix_y);

    diffMatrix_x = Lucee::Matrix<double>(nlocal, nlocal);
    diffMatrix_y = Lucee::Matrix<double>(nlocal, nlocal);

    for (unsigned i=0; i<nlocal; ++i)
      for (unsigned j=0; j<nlocal; ++j)
      {
        diffMatrix_x(i,j) = stiffMatrix_x(j,i);
        diffMatrix_y(i,j) = stiffMatrix_y(j,i);
      }

// NOTE: In the following calls, we need to keep getting the mass
// matrix again and again because the solve() method destroys
// it. Perhaps a better thing is to add new methods to the LinAlgebra
// file that computes the LU decomposition and then does the
// backsubstitution so all these inversions can be avoided.

    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

// multiply all matrices by M^-1
    nodalBasis->getMassMatrix(massMatrix);
    Lucee::solve(massMatrix, stiffMatrix_x);

    nodalBasis->getMassMatrix(massMatrix);
    Lucee::solve(massMatrix, stiffMatrix_y);

    nodalBasis->getMassMatrix(massMatrix);
    Lucee::solve(massMatrix, diffMatrix_x);

    nodalBasis->getMassMatrix(massMatrix);
    Lucee::solve(massMatrix, diffMatrix_y);
  }

  Lucee::UpdaterStatus 
  NodalPoissonBracketUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input/output arrays
    const Lucee::Field<2, double>& aCurr = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(1);

    Lucee::Field<2, double>& aNew = this->getOut<Lucee::Field<2, double> >(0);

// time-step
    double dt = t-this->getCurrTime();
// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// maximum CFL number used
    double cfla = 0.0;

// clear out contents of output field
    aNew = 0.0;

// compute contributions from volume integrals
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
      }
    }

// compute contributions for surface integrals. NOTE: In each
// direction, there is one extra edge than there are cells
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
      }
    }

    return Lucee::UpdaterStatus();
  }
  
  void
  NodalPoissonBracketUpdater::declareTypes()
  {
// takes two inputs (aOld, b)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// returns one output, (aNew)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}

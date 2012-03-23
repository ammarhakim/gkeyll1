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
#include <LcMathLib.h>
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

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around
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
      { // diff matrices are computed from transposed stiffness matrices
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
// get input arrays
    const Lucee::Field<2, double>& aCurr = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(1);
// get output arrays
    Lucee::Field<2, double>& aNew = this->getOut<Lucee::Field<2, double> >(0);

// time-step
    double dt = t-this->getCurrTime();
// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// maximum CFL number used
    double cfla = 0.0;

// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// space for potential and derivatives at nodes
    std::vector<double> phiK(nlocal), gradPhiK_x(nlocal), gradPhiK_y(nlocal);
    std::vector<double> flux_x(nlocal), flux_y(nlocal);
    double tv;

// various iterators
    Lucee::ConstFieldPtr<double> phiPtr = phi.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr = aCurr.createConstPtr();
    Lucee::FieldPtr<double> aNewPtr = aNew.createPtr();

// clear out contents of output field
    aNew = 0.0;

    double dx[2];

// compute contributions from volume integrals
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        
        nodalBasis->setIndex(ix, iy);
// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);

// compute gradient in X and Y directions
        calcGradient_x(phiK, gradPhiK_x);
        calcGradient_y(phiK, gradPhiK_y);

        aCurr.setPtr(aCurrPtr, ix, iy);
// compute x and y fluxes at each node
        for (unsigned k=0; k<nlocal; ++k)
        {
          flux_x[k] = gradPhiK_y[k]*aCurrPtr[k];
          flux_y[k] = -gradPhiK_x[k]*aCurrPtr[k];
        }

        aNew.setPtr(aNewPtr, ix, iy);
// multiply by stiffness matrix to give volume contribution
        for (unsigned i=0; i<nlocal; ++i)
        {
// NOTE: Can use BLAS to speed this up.
          tv = 0.0;
          for (unsigned j=0; j<nlocal; ++j)
            tv += (stiffMatrix_x(i,j)*flux_x[j] + stiffMatrix_y(i,j)*flux_y[j]);
          aNewPtr[i] += tv;
        }

// get grid spacing
        grid.setIndex(ix, iy);
        double dtdx = dt/grid.getDx(0);
        double dtdy = dt/grid.getDx(1);
// compute CFL number. THIS REALLY IS NOT CORRECT AS THE MAXIMUM SPEED
// SHOULD BE USED AND NOT JUST SPEED AT THE FIRST NODE. WILL FIX
// LATER, FOR NOW IT IS OKAY. Ammar, 3/22/2012
        cfla = Lucee::max3(cfla, dtdx*std::fabs(gradPhiK_y[0]), dtdy*std::fabs(gradPhiK_x[0]));
      }
    }

    if (cfla > cflm)
// time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

// compute contributions for surface integrals.
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
      }
    }

// perform final sweep, updating solution with forward Euler step
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        aNew.setPtr(aNewPtr, ix, iy);
        aCurr.setPtr(aCurrPtr, ix, iy);

        for (unsigned k=0; k<nlocal; ++k)
          aNewPtr[k] = aCurrPtr[k] + dt*aNewPtr[k];
      }
    }
    
    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
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

  void
  NodalPoissonBracketUpdater::calcGradient_x(const std::vector<double>& phiK, std::vector<double>& phiPrimeK)
  {
// NOTE: Can use BLAS to speed this up.

    unsigned n=phiK.size();
// gradient is multiplication by differentiation matrix
    for (unsigned i=0; i<n; ++i)
    {
      phiPrimeK[i] = 0.0;
      for (unsigned j=0; j<n; ++j)
        phiPrimeK[i] += diffMatrix_x(i,j)*phiK[j];
    }
  }

  void
  NodalPoissonBracketUpdater::calcGradient_y(const std::vector<double>& phiK, std::vector<double>& phiPrimeK)
  {
// NOTE: Can use BLAS to speed this up.

    unsigned n=phiK.size();
// gradient is multiplication by differentiation matrix
    for (unsigned i=0; i<n; ++i)
    {
      phiPrimeK[i] = 0.0;
      for (unsigned j=0; j<n; ++j)
        phiPrimeK[i] += diffMatrix_y(i,j)*phiK[j];
    }
  }
}

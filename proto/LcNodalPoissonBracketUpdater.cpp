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
  static const unsigned UPWIND = 0;
  static const unsigned CENTRAL = 1;

// set id for module system
  const char *NodalPoissonBracketUpdater::id = "PoissonBracket";

  NodalPoissonBracketUpdater::NodalPoissonBracketUpdater()
    : UpdaterIfc()
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

    fluxType = UPWIND;
    if (tbl.hasString("fluxType"))
    {
      if (tbl.getString("fluxType") == "upwind")
        fluxType = UPWIND;
      else if (tbl.getString("fluxType") == "central")
        fluxType = CENTRAL;
      else
      {
        Lucee::Except lce("NodalPoissonBracketUpdater::readInput: 'fluxType' ");
        lce << tbl.getString("fluxType") << " is not valid";
        throw lce;
      }
    }
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

    unsigned nlocal = nodalBasis->getNumNodes();

// get node numbers on each lower and upper edges
    for (unsigned dir=0; dir<2; ++dir)
    {
      lowerNodeNums[dir].nums.resize(nodalBasis->getNumSurfLowerNodes(dir));
      nodalBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

      upperNodeNums[dir].nums.resize(nodalBasis->getNumSurfUpperNodes(dir));
      nodalBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);
    }

// space for mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<2; ++dir)
    {
// get stiffness matrice
      stiffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, stiffMatrix[dir].m);

// calculate differentiation matrix
      diffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      for (unsigned i=0; i<nlocal; ++i)
        for (unsigned j=0; j<nlocal; ++j)
// diff matrices are computed from transposed stiffness matrices
          diffMatrix[dir].m(i,j) = stiffMatrix[dir].m(j,i);

// multiply matrices by inverse of mass matrix
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, stiffMatrix[dir].m);

      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, diffMatrix[dir].m);

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
// number of nodes on each face
    unsigned nFace = 2;  // NOTE: THIS NEEDS TO BE GOTTEN FROM nodalBasis

// space for potential and derivatives at nodes
    std::vector<double> phiK(nlocal), gradPhiK_x(nlocal), gradPhiK_y(nlocal);
    std::vector<double> flux_x(nlocal), flux_y(nlocal);
    std::vector<double> fdotn(nFace);
    double tv;

// various iterators
    Lucee::ConstFieldPtr<double> phiPtr = phi.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_l = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_r = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_t = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_b = aCurr.createConstPtr();
    Lucee::FieldPtr<double> aNewPtr = aNew.createPtr();

// clear out contents of output field
    aNew = 0.0;

    double dx[2];

// compute contributions from volume integrals
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        aCurr.setPtr(aCurrPtr, ix, iy);
        aCurr.setPtr(aCurrPtr_l, ix-1, iy); // left cell
        aCurr.setPtr(aCurrPtr_r, ix+1, iy); // right cell
        aCurr.setPtr(aCurrPtr_t, ix, iy+1); // top cell
        aCurr.setPtr(aCurrPtr_b, ix, iy-1); // bottom cell

        aNew.setPtr(aNewPtr, ix, iy);
        nodalBasis->setIndex(ix, iy);

// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);

// compute gradient in X and Y directions
        calcGradient_x(phiK, gradPhiK_x);
        calcGradient_y(phiK, gradPhiK_y);

// compute x and y fluxes at each node
        for (unsigned k=0; k<nlocal; ++k)
        {
          flux_x[k] = gradPhiK_y[k]*aCurrPtr[k];
          flux_y[k] = -gradPhiK_x[k]*aCurrPtr[k];
        }

// multiply by stiffness matrix to give volume contribution
        matVec(1.0, stiffMatrix[0].m, flux_x, 1.0, &aNewPtr[0]);
        matVec(1.0, stiffMatrix[1].m, flux_y, 1.0, &aNewPtr[0]);

// compute contribution from edge integrals

// edge X-lower has nodes (1,4)
        fdotn[0] = -getUpwindFlux(gradPhiK_y[0], aCurrPtr_l[1], aCurrPtr[0]);
        fdotn[1] = -getUpwindFlux(gradPhiK_y[3], aCurrPtr_l[2], aCurrPtr[3]);

// multiply by lifting matrix
        matVec(-1.0, lowerLift[0].m, fdotn, 1.0, &aNewPtr[0]);

// edge X-upper has nodes (2,3)
        fdotn[0] = getUpwindFlux(gradPhiK_y[1], aCurrPtr[1], aCurrPtr_r[0]);
        fdotn[1] = getUpwindFlux(gradPhiK_y[2], aCurrPtr[2], aCurrPtr_r[3]);

// multiply by lifting matrix
        matVec(-1.0, upperLift[0].m, fdotn, 1.0, &aNewPtr[0]);

// edge Y-lower has nodes (1,2)
        fdotn[0] = -getUpwindFlux(-gradPhiK_x[0], aCurrPtr_b[3], aCurrPtr[0]);
        fdotn[1] = -getUpwindFlux(-gradPhiK_x[1], aCurrPtr_b[2], aCurrPtr[1]);

// multiply by lifting matrix
        matVec(-1.0, lowerLift[1].m, fdotn, 1.0, &aNewPtr[0]);

// edge Y-upper has nodes (4,3)
        fdotn[0] = getUpwindFlux(-gradPhiK_x[3], aCurrPtr[3], aCurrPtr_t[0]);
        fdotn[1] = getUpwindFlux(-gradPhiK_x[2], aCurrPtr[2], aCurrPtr_t[1]);

// multiply by lifting matrix
        matVec(-1.0, upperLift[1].m, fdotn, 1.0, &aNewPtr[0]);

// get grid spacing
        grid.setIndex(ix, iy);
        double dtdx = dt/grid.getDx(0);
        double dtdy = dt/grid.getDx(1);
// compute CFL number.
        for (unsigned n=0; n<nlocal; ++n)
          cfla = Lucee::max3(cfla, dtdx*std::fabs(gradPhiK_y[n]), dtdy*std::fabs(gradPhiK_x[n]));
      }
    }

    if (cfla > cflm)
// time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

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
  NodalPoissonBracketUpdater::calcGradient_x(std::vector<double>& phiK, std::vector<double>& phiPrimeK)
  {
    matVec(1.0, diffMatrix[0].m, phiK, 0.0, &phiPrimeK[0]);
  }

  void
  NodalPoissonBracketUpdater::calcGradient_y(std::vector<double>& phiK, std::vector<double>& phiPrimeK)
  {
    matVec(1.0, diffMatrix[1].m, phiK, 0.0, &phiPrimeK[0]);
  }

  double
  NodalPoissonBracketUpdater::getUpwindFlux(double u, double chil, double chir)
  {
    if (fluxType == UPWIND)
    {
      if (u > 0.0)
        return u*chil;
      else
        return u*chir;
    }
    else if (fluxType == CENTRAL)
      return u*0.5*(chil+chir);

    return 0.0;
  }

  void 
  NodalPoissonBracketUpdater::matVec(double m, const Lucee::Matrix<double>& mat,
    const std::vector<double>& vec, double v, double *out)
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
}

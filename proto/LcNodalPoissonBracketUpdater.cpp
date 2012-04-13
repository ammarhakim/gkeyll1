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
// number of nodes on each face (WARNING: assumption here is that all
// faces have same number of nodes)
    unsigned nFace = nodalBasis->getNumSurfLowerNodes(0);

// space for various quantities
    std::vector<double> phiK(nlocal), flux(nlocal);
    std::vector<double> fdotn(nFace);

    NodalPoissonBracketUpdater::NodeSpeed speeds[2];
    for (unsigned dir=0; dir<2; ++dir)
      speeds[dir].s.resize(nlocal);

// various iterators
    Lucee::ConstFieldPtr<double> phiPtr = phi.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_n = aCurr.createConstPtr();
    Lucee::FieldPtr<double> aNewPtr = aNew.createPtr();

// clear out contents of output field
    aNew = 0.0;

    double dx[2];
    int idx[2];

// compute contributions from volume integrals
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        aCurr.setPtr(aCurrPtr, ix, iy);
        aNew.setPtr(aNewPtr, ix, iy);
        nodalBasis->setIndex(ix, iy);

// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);

// compute speeds
        calcSpeeds(phiK, speeds);

// compute fluxes at each interior node and accumulate contribution to
// volume integral
        for (unsigned dir=0; dir<2; ++dir)
        {
          for (unsigned k=0; k<nlocal; ++k)
            flux[k] = speeds[dir].s[k]*aCurrPtr[k];
          matVec(1.0, stiffMatrix[dir].m, flux, 1.0, &aNewPtr[0]);
        }

// compute contribution from edge integrals on lower edges
        for (unsigned dir=0; dir<2; ++dir)
        {
          idx[0] = ix; idx[1] = iy;
          idx[dir] += -1; // lower side

// set pointer to cell on lower side in this direction
          aCurr.setPtr(aCurrPtr_n, idx[0], idx[1]);

// contribution from lower edge in direction
          for (unsigned f=0; f<nFace; ++f)
          {
            unsigned fn = lowerNodeNums[dir].nums[f]; // node number in cell on right
            unsigned fn_n = upperNodeNums[dir].nums[f]; // node number in cell on left

// compute upwind flux (lower edges contribution is -ve)
            fdotn[f] = -getUpwindFlux(speeds[dir].s[fn],
              aCurrPtr_n[fn_n], aCurrPtr[fn]);
          }
// multiply by lifting matrix to compute contribution
          matVec(-1.0, lowerLift[dir].m, fdotn, 1.0, &aNewPtr[0]);
        }

// compute contribution from edge integrals on upper edges
        for (unsigned dir=0; dir<2; ++dir)
        {
          idx[0] = ix; idx[1] = iy;
          idx[dir] += 1; // upper side

// set pointer to cell on lower side in this direction
          aCurr.setPtr(aCurrPtr_n, idx[0], idx[1]);

// contribution from lower edge in direction
          for (unsigned f=0; f<nFace; ++f)
          {
            unsigned fn = upperNodeNums[dir].nums[f]; // node number in cell on left
            unsigned fn_n = lowerNodeNums[dir].nums[f]; // node number in cell on right

// compute upwind flux (upper edges contribution is +ve)
            fdotn[f] = getUpwindFlux(speeds[dir].s[fn],
              aCurrPtr[fn], aCurrPtr_n[fn_n]);
          }
// multiply by lifting matrix to compute contribution
          matVec(-1.0, upperLift[dir].m, fdotn, 1.0, &aNewPtr[0]);
        }

// get grid spacing
        grid.setIndex(ix, iy);
        double dtdx = dt/grid.getDx(0);
        double dtdy = dt/grid.getDx(1);
// compute CFL number.
        for (unsigned n=0; n<nlocal; ++n)
        {
          cfla = Lucee::max3(cfla, dtdx*std::fabs(speeds[0].s[n]), 
            dtdy*std::fabs(speeds[1].s[n]));
        }
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

  void
  NodalPoissonBracketUpdater::calcSpeeds(std::vector<double>& phiK, NodeSpeed speeds[2])
  {
// ux = d phi / dy
    matVec(1.0, diffMatrix[1].m, phiK, 0.0, &speeds[0].s[0]);
// uy = - d phi / dx
    matVec(-1.0, diffMatrix[0].m, phiK, 0.0, &speeds[1].s[0]);
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

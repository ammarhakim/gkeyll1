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

// If any value is negative, sets it to zero  
  static
  void applyPositivityFix(std::vector<double>& vals)
  {
    for (unsigned k=0; k<vals.size(); ++k)
      if (vals[k] < 0.0) vals[k] = 0.0;
  }

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

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("NodalPoissonBracketUpdater::readInput: Must specify element to use using 'basis'");

    cfl = tbl.getNumber("cfl");
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

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    hamilNodesShared = true;
    if (tbl.hasBool("hamilNodesShared"))
      hamilNodesShared = tbl.getBool("hamilNodesShared");

    applyPositivityFluxFix = false;
    if (tbl.hasBool("applyPositivityFix"))
      applyPositivityFluxFix = tbl.getBool("applyPositivityFix");

// zeroFluxOffset allows implementation of zero-flux BCs in specified
// directions. This is needed to avoid "leakage" of particles when
// doing Vlasov solves, as velocity space direction is bounded.
     
    lowerZeroFluxOffset[0] = 0; lowerZeroFluxOffset[1] = 0; // by default no offsets
    upperZeroFluxOffset[0] = 0; upperZeroFluxOffset[1] = 0; // by default no offsets

    std::vector<double> zfd; // zero-flux directions
    if (tbl.hasNumVec("zeroFluxDirections"))
      zfd = tbl.getNumVec("zeroFluxDirections");
    for (unsigned i=0; i<zfd.size(); ++i)
    {
      if ((int)zfd[i]>1)
        throw Lucee::Except("NodalPoissonBracketUpdater::readInput: Directions in 'zeroFluxDirections' should be 0 or 1");
      lowerZeroFluxOffset[(int)zfd[i]] = 1;
      upperZeroFluxOffset[(int)zfd[i]] = 1;
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

    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<2; ++dir)
    {
// get stiffness matrix
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
    }

    unsigned nVolQuad = nodalBasis->getNumGaussNodes();
    volQuad.reset(nVolQuad, nlocal);

// get data needed for Gaussian quadrature
    nodalBasis->getGaussQuadData(volQuad.interpMat.m, volQuad.ords.m, volQuad.weights);

    for (unsigned dir=0; dir<2; ++dir)
    {
// compute differentiation matrices that compute derivatives at
// quadrature nodes
      Lucee::accumulate(volQuad.pDiffMatrix[dir].m, volQuad.interpMat.m, diffMatrix[dir].m);

      mGradPhi[dir].m = Lucee::Matrix<double>(nlocal, nVolQuad);
// compute gradient of basis functions at quadrature nodes (this is
// just differentiation matrix at quadrature nodes)
      for (unsigned i=0; i<nlocal; ++i)
        for (unsigned j=0; j<nVolQuad; ++j)
// diff matrices are computed from transposed stiffness matrices
          mGradPhi[dir].m(i,j) = volQuad.pDiffMatrix[dir].m(j,i);

// multiply by inverse mass matrix at this point
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, mGradPhi[dir].m);
    }

    unsigned nSurfQuad = nodalBasis->getNumSurfGaussNodes();
// get data for surface quadrature
    for (unsigned dir=0; dir<2; ++dir)
    {
      surfLowerQuad[dir].reset(nSurfQuad, nlocal);
// lower surface data
      nodalBasis->getSurfLowerGaussQuadData(dir, surfLowerQuad[dir].interpMat.m,
        surfLowerQuad[dir].ords.m, surfLowerQuad[dir].weights);

      surfUpperQuad[dir].reset(nSurfQuad, nlocal);
// upper surface data
      nodalBasis->getSurfUpperGaussQuadData(dir, surfUpperQuad[dir].interpMat.m,
        surfUpperQuad[dir].ords.m, surfUpperQuad[dir].weights);
    }

    for (unsigned dir=0; dir<2; ++dir)
    { // dir is direction normal to face

// compute differentiation matrices that compute derivatives at
// surface quadrature nodes
      for (unsigned d=0; d<2; ++d)
      { // d direction of derivative
        Lucee::accumulate(surfLowerQuad[dir].pDiffMatrix[d].m, 
          surfLowerQuad[dir].interpMat.m, diffMatrix[d].m);

        Lucee::accumulate(surfUpperQuad[dir].pDiffMatrix[d].m, 
          surfUpperQuad[dir].interpMat.m, diffMatrix[d].m);
      }

      mSurfLowerPhi[dir].m = Lucee::Matrix<double>(nlocal, nSurfQuad);
// compute basis functions at surface quadrature nodes
      for (unsigned i=0; i<nlocal; ++i)
        for (unsigned j=0; j<nSurfQuad; ++j)
          mSurfLowerPhi[dir].m(i,j) = surfLowerQuad[dir].interpMat.m(j,i);

// multiply by inverse mass matrix at this point
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, mSurfLowerPhi[dir].m);

      mSurfUpperPhi[dir].m = Lucee::Matrix<double>(nlocal, nSurfQuad);
// compute basis functions at surface quadrature nodes
      for (unsigned i=0; i<nlocal; ++i)
        for (unsigned j=0; j<nSurfQuad; ++j)
          mSurfUpperPhi[dir].m(i,j) = surfUpperQuad[dir].interpMat.m(j,i);

// multiply by inverse mass matrix at this point
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, mSurfUpperPhi[dir].m);
    }

// ensure that zero-flux BCs are applied only if local rank owns the
// skin cell
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();
    for (unsigned d=0; d<2; ++d)
    {
      if (localRgn.getLower(d) != globalRgn.getLower(d))
        lowerZeroFluxOffset[d] = 0; // not owned by us, so ignore
      if (localRgn.getUpper(d) != globalRgn.getUpper(d))
        upperZeroFluxOffset[d] = 0; // not owned by us, so ignore
    }      
  }

  Lucee::UpdaterStatus 
  NodalPoissonBracketUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    const Lucee::Field<2, double>& aCurr = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(1);

    Lucee::Field<2, double>& aNew = this->getOut<Lucee::Field<2, double> >(0);

    double dt = t-this->getCurrTime();
// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    double cfla = 0.0; // maximum CFL number used

    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned nVolQuad = nodalBasis->getNumGaussNodes();
    unsigned nFace = nodalBasis->getNumSurfLowerNodes(0);
    unsigned nSurfQuad = nodalBasis->getNumSurfGaussNodes();

    std::vector<double> phiK(nlocal);
    std::vector<double> chiQuad_l(nSurfQuad), chiQuad_r(nSurfQuad), chiUpwind(nSurfQuad);
    std::vector<double> udotn(nSurfQuad), uflux(nSurfQuad), fdotn(nSurfQuad);
    std::vector<double> quadChi(nVolQuad);

    NodalPoissonBracketUpdater::NodeSpeed speeds[2], quadSpeeds[2];
    for (unsigned dir=0; dir<2; ++dir)
    {
      speeds[dir].s.resize(nVolQuad);
      quadSpeeds[dir].s.resize(nVolQuad);
    }

    Lucee::ConstFieldPtr<double> phiPtr = phi.createConstPtr();
    Lucee::ConstFieldPtr<double> phiPtr_l = phi.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_l = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_n = aCurr.createConstPtr();
    Lucee::FieldPtr<double> aNewPtr = aNew.createPtr();
    Lucee::FieldPtr<double> aNewPtr_l = aNew.createPtr();

    aNew = 0.0; // use aNew to store increment initially

    double dx[2];
    int idx[2];

// compute contributions from volume integrals
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        aCurr.setPtr(aCurrPtr, ix, iy);
        aNew.setPtr(aNewPtr, ix, iy);

// extract potential at this location
        if (hamilNodesShared)
        {
          nodalBasis->setIndex(ix, iy);        
          nodalBasis->extractFromField(phi, phiK);
        }
        else
        {
          phi.setPtr(phiPtr, ix, iy);
          for (unsigned k=0; k<nlocal; ++k) phiK[k] = phiPtr[k];
        }

// compute speeds
        calcSpeedsAtQuad(phiK, quadSpeeds);

// interpolate vorticity to quadrature points
        matVec(1.0, volQuad.interpMat.m, &aCurrPtr[0], 0.0, &quadChi[0]);
// zero out negative values if positivity fix for flux is needed
        if (applyPositivityFluxFix)
          applyPositivityFix(quadChi);

// compute fluxes at each interior node and accumulate contribution to
// volume integral
        for (unsigned dir=0; dir<2; ++dir)
        {
          for (unsigned k=0; k<nlocal; ++k)
          {
// loop over quadrature points, accumulating contribution
            for (unsigned qp=0; qp<nVolQuad; ++qp)
              aNewPtr[k] += volQuad.weights[qp]*mGradPhi[dir].m(k,qp)*quadSpeeds[dir].s[qp]*quadChi[qp];
          }
        }

        grid.setIndex(ix, iy);
        double dtdx = dt/grid.getDx(0);
        double dtdy = dt/grid.getDx(1);
// compute CFL number.
        for (unsigned n=0; n<nVolQuad; ++n)
        {
          cfla = Lucee::max3(cfla, dtdx*std::fabs(quadSpeeds[0].s[n]),
            dtdy*std::fabs(quadSpeeds[1].s[n]));
        }
      }
    }

    if (cfla > cflm)
// time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

// these strange arrays hold the direction and signs to compute
// gradients of phi given the orientation of an edge. In X-direction
// we need the Y gradients of phi, while in Y-direction we need to
// compute the X gradients of phi.
    unsigned gradDir[2] = {1, 0};
    int signs[2] = {1.0, -1.0};

// perform surface integrals
    for (unsigned dir=0; dir<2; ++dir)
    {
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<2> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice. (We need to make sure that flux
// is computed for one edge outside the domain interior, accounting
// for the fact that we may not want to compute fluxes from the
// outermost edges [zero-flux BCs])
      int sliceLower = localRgn.getLower(dir)+lowerZeroFluxOffset[dir];
      int sliceUpper = localRgn.getUpper(dir)+1-upperZeroFluxOffset[dir];

      int idx[2], idxl[2];
      int ix, iy;
// loop over each 1D slice
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxl);

// loop over each edge in slice
        for (int i=sliceLower; i<sliceUpper; ++i)
        {
          idx[dir] = i; // cell right of edge
          idxl[dir] = i-1; // cell left of edge

          aCurr.setPtr(aCurrPtr, idx);
          aCurr.setPtr(aCurrPtr_l, idxl);

// NOTE: We extract potential from the cell left of the edge. We can
// also use the cell right of the edge as u.nhat is
// continous. However, the right cell on the upper domain boundary
// does not have complete data due to the strange way in which the
// nodal data is stored for continous FE fields. This is a very subtle
// issue and eventually needs to be fixed. (This is not a problem if
// the Hamiltonian is in a DG field (even if it is continuous), which
// should be the prefered way to do things now. AH: 8/25/2015)
          if (hamilNodesShared)
          {
            nodalBasis->setIndex(idxl);
            nodalBasis->extractFromField(phi, phiK);
          }
          else
          {
            phi.setPtr(phiPtr, idxl);
            for (unsigned k=0; k<nlocal; ++k) phiK[k] = phiPtr[k];
          }

// now compute u.nhat on this edge (it is the upper edge of the cell
// left of the edge)
          matVec(signs[dir], surfUpperQuad[dir].pDiffMatrix[gradDir[dir]].m, &phiK[0],
            0.0, &udotn[0]);

// interpolate vorticity to the quadrature nodes (left side is upper
// edge of left cell, while right side is lower edge of right cell)
          matVec(1.0, surfUpperQuad[dir].interpMat.m, &aCurrPtr_l[0], 0.0, &chiQuad_l[0]);
          matVec(1.0, surfLowerQuad[dir].interpMat.m, &aCurrPtr[0], 0.0, &chiQuad_r[0]);

// zero out negative values if positivity fix for flux is needed
          if (applyPositivityFluxFix)
          {
            applyPositivityFix(chiQuad_l);
            applyPositivityFix(chiQuad_r);
          }

          for (unsigned qp=0; qp<nSurfQuad; ++qp)
            uflux[qp] = getUpwindFlux(udotn[qp], chiQuad_l[qp], chiQuad_r[qp]);

// at this point we have the flux at the edge. We need to accumulate
// its contribution to the cells connected to the edge
          aNew.setPtr(aNewPtr, idx);
          aNew.setPtr(aNewPtr_l, idxl);

// perform surface integration
          for (unsigned k=0; k<nlocal; ++k)
          {
            for (unsigned qp=0; qp<nSurfQuad; ++qp)
              aNewPtr[k] += surfLowerQuad[dir].weights[qp]*mSurfLowerPhi[dir].m(k,qp)*uflux[qp];
          }

// perform surface integration
          for (unsigned k=0; k<nlocal; ++k)
          {
            for (unsigned qp=0; qp<nSurfQuad; ++qp)
              aNewPtr_l[k] += -surfUpperQuad[dir].weights[qp]*mSurfUpperPhi[dir].m(k,qp)*uflux[qp];
          }
        }
      }
    }

// NOTE: If only calculation of increments are requested, the final
// Euler update is not performed. This means that the multiplication
// of the DG RHS with dt is not done, something to keep in mind if
// using the increment in time-dependent update.

    if (onlyIncrement == false)
    { // only do this if full update is requested

// Perform final sweep, updating solution with forward Euler step
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
  NodalPoissonBracketUpdater::calcSpeeds(std::vector<double>& phiK, NodeSpeed speeds[2])
  {
// ux = d phi / dy
    matVec(1.0, diffMatrix[1].m, &phiK[0], 0.0, &speeds[0].s[0]);
// uy = - d phi / dx
    matVec(-1.0, diffMatrix[0].m, &phiK[0], 0.0, &speeds[1].s[0]);
  }

  void
  NodalPoissonBracketUpdater::calcSpeedsAtQuad(std::vector<double>& phiK, NodeSpeed speeds[2])
  {
// ux = d phi / dy
    matVec(1.0, volQuad.pDiffMatrix[1].m, &phiK[0], 0.0, &speeds[0].s[0]);
// uy = - d phi / dx
    matVec(-1.0, volQuad.pDiffMatrix[0].m, &phiK[0], 0.0, &speeds[1].s[0]);
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
}

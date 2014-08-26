/**
 * @file	LcPoissonBracketUpdater.cpp
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
#include <LcPoissonBracketUpdater.h>

namespace Lucee
{
  static const unsigned UPWIND = 0;
  static const unsigned CENTRAL = 1;

// set id for module system
  template <> const char *PoissonBracketUpdater<1>::id = "PoissonBracket1D";
  template <> const char *PoissonBracketUpdater<2>::id = "PoissonBracket2D";
  template <> const char *PoissonBracketUpdater<3>::id = "PoissonBracket3D";
  template <> const char *PoissonBracketUpdater<4>::id = "PoissonBracket4D";
  template <> const char *PoissonBracketUpdater<5>::id = "PoissonBracket5D";

  template <unsigned NDIM>
  PoissonBracketUpdater<NDIM>::PoissonBracketUpdater()
    : UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void 
  PoissonBracketUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("PoissonBracketUpdater::readInput: Must specify element to use using 'basis'");

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
        Lucee::Except lce("PoissonBracketUpdater::readInput: 'fluxType' ");
        lce << tbl.getString("fluxType") << " is not valid";
        throw lce;
      }
    }

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }

  template <unsigned NDIM>
  void 
  PoissonBracketUpdater<NDIM>::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    
    unsigned nlocal = nodalBasis->getNumNodes();

    // set index to first location in grid (this is okay as in this
    // updater we are assuming grid is uniform)
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step();
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

    // Store mass matrix inverse
    Lucee::Matrix<double> tempMass(nlocal, nlocal);
    nodalBasis->getMassMatrix(tempMass);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMass, massMatrix);
    massMatrixInv = massMatrix.inverse();

    // Store grad stiffness matrix in each direction
    gradStiffnessMatrix.resize(NDIM);
    for (int dir = 0; dir < NDIM; dir++)
    {
      Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, tempMatrix);
      gradStiffnessMatrix[dir] = Eigen::MatrixXd(nlocal, nlocal);

      copyLuceeToEigen(tempMatrix, gradStiffnessMatrix[dir]);
    }

    // get node numbers on each lower and upper edges
    lowerNodeNums.resize(NDIM);
    upperNodeNums.resize(NDIM);
    for (int dir = 0; dir < NDIM; dir++)
    {
      lowerNodeNums[dir].resize(nodalBasis->getNumSurfLowerNodes(dir));
      nodalBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir]);

      upperNodeNums[dir].resize(nodalBasis->getNumSurfUpperNodes(dir));
      nodalBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir]);
    }

    // get data needed for Gaussian quadrature
    int nVolQuad = nodalBasis->getNumGaussNodes();
    volQuad.reset(nVolQuad, nlocal);
    std::vector<double> volWeights(nVolQuad);

    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NDIM);
    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex];
    
    std::cout << "weights" << std::endl << volQuad.weights << std::endl;

    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
    copyLuceeToEigen(tempVolCoords, volQuad.coordMat);

    for (int dir = 0; dir < NDIM; dir++)
    {
      // Each row is a quadrature point; each column is a basis function with derivative applied
      Eigen::MatrixXd derivMatrix = volQuad.interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

      // Copy derivatives to different kind type of matrix
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        volQuad.pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);
    }

    /*
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
    */
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  PoissonBracketUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& aCurr = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& phi = this->getInp<Lucee::Field<NDIM, double> >(1);

    Lucee::Field<NDIM, double>& aNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    double dt = t-this->getCurrTime();
// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    double cfla = 0.0; // maximum CFL number used

    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned nVolQuad = nodalBasis->getNumGaussNodes();
    unsigned nFace = nodalBasis->getNumSurfLowerNodes(0);
    unsigned nSurfQuad = nodalBasis->getNumSurfGaussNodes();

    std::vector<double> phiK(nlocal);
    std::vector<double> chiQuad_l(nSurfQuad), chiQuad_r(nSurfQuad), chiUpwind(nSurfQuad);
    std::vector<double> udotn(nSurfQuad), uflux(nSurfQuad), fdotn(nSurfQuad);
    std::vector<double> quadChi(nVolQuad);

    NodeSpeed speeds[2], quadSpeeds[2];
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

    double dx[NDIM];

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    // contributions from volume integrals
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      aCurr.setPtr(aCurrPtr, idx);
      aNew.setPtr(aNewPtr, idx);

      Eigen::MatrixXd alpha(NDIM, nVolQuad);

      // Get a vector of f at quad points
      Eigen::VectorXd fVec(nlocal);
      for (int i = 0; i < nlocal; i++)
        fVec(i) = aCurrPtr[i];
      Eigen::VectorXd fAtQuad = volQuad.interpMat*fVec;

      // Coefficient-wise multiply each row of alpha by fAtQuad
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        for (int quadIndex = 0; quadIndex < fAtQuad.size(); quadIndex++)
          alpha(dimIndex, quadIndex) *= fAtQuad(quadIndex);

      Eigen::VectorXd resultVector(nlocal);
      // Compute grad of a test function, multiply it with alpha and f, then multiply
      // with weights and sum to compute integral.
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        resultVector(basisIndex) = (volQuad.pDiffMatrix[basisIndex].cwiseProduct(alpha)*volQuad.weights).sum();

      // Multiply resultVector with massMatrixInv to calculate aNew updates
    }

    // TEST CODE

    // contributions from surface integrals
    for (int dir = 0; dir < NDIM; dir++)
    {
      // create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));
      // lower and upper bounds of 1D slice. (We need to make sure that flux
      // is computed for one edge outside domain interior)
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

      int idxr[NDIM];
      int idxl[NDIM];

      // loop over each 1D slice
      while (seq.step())
      {
        seq.fillWithIndex(idxr);
        seq.fillWithIndex(idxl);
        // loop over each edge
        for (int i = sliceLower; i < sliceUpper; i++)
        { 
          idxr[dir] = i;
          idxl[dir] = i-1;

          grid.setIndex(idxl);
          double dxL = grid.getDx(dir);
          grid.setIndex(idxr);
          double dxR = grid.getDx(dir);
        }
      }
    }
/*
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
        calcSpeedsAtQuad(phiK, quadSpeeds);

// interpolate vorticity to quadrature points
        matVec(1.0, volQuad.interpMat.m, &aCurrPtr[0], 0.0, &quadChi[0]);

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
        for (unsigned n=0; n<nlocal; ++n)
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

// lower and upper bounds of 1D slice. (We need to make sure that the
// flux is computed for one edge outside the domain interior)
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

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
// issue and eventually needs to be fixed.
          nodalBasis->setIndex(idxl);
          nodalBasis->extractFromField(phi, phiK);

// now compute u.nhat on this edge (it is the upper edge of the cell
// left of the edge)
          matVec(signs[dir], surfUpperQuad[dir].pDiffMatrix[gradDir[dir]].m, &phiK[0],
            0.0, &udotn[0]);

// interpolate vorticity to the quadrature nodes (left side is upper
// edge of left cell, while right side is lower edge of right cell)
          matVec(1.0, surfUpperQuad[dir].interpMat.m, &aCurrPtr_l[0], 0.0, &chiQuad_l[0]);
          matVec(1.0, surfLowerQuad[dir].interpMat.m, &aCurrPtr[0], 0.0, &chiQuad_r[0]);

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
    */
    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }
  
  template <unsigned NDIM>
  void
  PoissonBracketUpdater<NDIM>::declareTypes()
  {
// takes two inputs (aOld, b)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// returns one output, (aNew)
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  PoissonBracketUpdater<NDIM>::calcSpeeds(std::vector<double>& phiK, NodeSpeed speeds[2])
  {
// ux = d phi / dy
    //matVec(1.0, diffMatrix[1], &phiK[0], 0.0, &speeds[0].s[0]);
// uy = - d phi / dx
    //matVec(-1.0, diffMatrix[0], &phiK[0], 0.0, &speeds[1].s[0]);
  }

  template <unsigned NDIM>
  void
  PoissonBracketUpdater<NDIM>::calcSpeedsAtQuad(std::vector<double>& phiK, NodeSpeed speeds[2])
  {
// ux = d phi / dy
    //matVec(1.0, volQuad.pDiffMatrix[1].m, &phiK[0], 0.0, &speeds[0].s[0]);
// uy = - d phi / dx
    //matVec(-1.0, volQuad.pDiffMatrix[0].m, &phiK[0], 0.0, &speeds[1].s[0]);
  }

  template <unsigned NDIM>
  double
  PoissonBracketUpdater<NDIM>::getUpwindFlux(double u, double chil, double chir)
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

  template <unsigned NDIM>
  void 
  PoissonBracketUpdater<NDIM>::matVec(double m, const Lucee::Matrix<double>& mat,
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

  template <unsigned NDIM>
  void
  PoissonBracketUpdater<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class PoissonBracketUpdater<1>;
  template class PoissonBracketUpdater<2>;
  template class PoissonBracketUpdater<3>;
  template class PoissonBracketUpdater<4>;
  template class PoissonBracketUpdater<5>;
}

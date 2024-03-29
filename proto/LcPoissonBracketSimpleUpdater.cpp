/**
 * @file	LcPoissonBracketSimpleUpdater.cpp
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 * This version assumes there is no spatially-varying transformation jacobian, allowing
 * for some speedups (?).
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMathLib.h>
#include <LcPoissonBracketSimpleUpdater.h>

namespace Lucee
{
  static const unsigned UPWIND = 0;
  static const unsigned CENTRAL = 1;
  static const unsigned DOWNWIND = 2;

// set id for module system
  template <> const char *PoissonBracketSimpleUpdater<1>::id = "PoissonBracketSimple1D";
  template <> const char *PoissonBracketSimpleUpdater<2>::id = "PoissonBracketSimple2D";
  template <> const char *PoissonBracketSimpleUpdater<3>::id = "PoissonBracketSimple3D";
  template <> const char *PoissonBracketSimpleUpdater<4>::id = "PoissonBracketSimple4D";
  template <> const char *PoissonBracketSimpleUpdater<5>::id = "PoissonBracketSimple5D";

  template <unsigned NDIM>
  PoissonBracketSimpleUpdater<NDIM>::PoissonBracketSimpleUpdater()
    : UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void 
  PoissonBracketSimpleUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("PoissonBracketSimpleUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasObject<Lucee::PoissonBracketEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::PoissonBracketEquation>("equation");
    else
    {
      Lucee::Except lce("PoissonBracketSimpleUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }

    cfl = tbl.getNumber("cfl");
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around

    fluxType = UPWIND;
    if (tbl.hasString("fluxType"))
    {
      if (tbl.getString("fluxType") == "upwind")
        fluxType = UPWIND;
      else if (tbl.getString("fluxType") == "central")
        fluxType = CENTRAL;
      else if (tbl.getString("fluxType") == "downwind")
        fluxType = DOWNWIND;
      else
      {
        Lucee::Except lce("PoissonBracketSimpleUpdater::readInput: 'fluxType' ");
        lce << tbl.getString("fluxType") << " is not valid";
        throw lce;
      }
    }

    jacobianFactor = 1.0;
    if (tbl.hasNumber("jacobianFactor"))
      jacobianFactor = tbl.getNumber("jacobianFactor");

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    // directions to update
    if (tbl.hasNumVec("updateDirections"))
    {
      std::vector<double> tempDirs = tbl.getNumVec("updateDirections");
      for (int i = 0; i < std::min<int>(NDIM, tempDirs.size()); i++)
      {
        int d = (int) tempDirs[i];
        if (d < NDIM)
          updateDirs.push_back(d);
        else
          throw Lucee::Except("PoissonBracketSimpleUpdater::readInput: updateDirections must be a table with each element < NDIM");
      }
    }
    else
    {
      for (int i = 0; i < NDIM; i++)
        updateDirs.push_back(i);
    }

    zeroFluxFlags = std::vector<bool>(NDIM);
    // Default is no zero flux bcs in any directions
    for (int i = 0; i < NDIM; i++)
      zeroFluxFlags[i] = false;

    if (tbl.hasNumVec("zeroFluxDirections"))
    {
      std::vector<double> tempDirs = tbl.getNumVec("zeroFluxDirections");
      for (int i = 0; i < std::min<int>(NDIM, tempDirs.size()); i++)
      {
        int d = (int) tempDirs[i];

        if (d < NDIM)
          zeroFluxFlags[d] = true;
        else
          throw Lucee::Except("PoissonBracketSimpleUpdater::readInput: zeroFluxDirections must be a table with each element < NDIM");
      }
    }
  }

  template <unsigned NDIM>
  void 
  PoissonBracketSimpleUpdater<NDIM>::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    
    int nlocal = nodalBasis->getNumNodes();

    // set index to first location in grid (this is okay as in this
    // updater we are assuming grid is uniform)
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step();
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    grid.setIndex(idx);

    // Store mass matrix inverse
    Lucee::Matrix<double> tempMass(nlocal, nlocal);
    nodalBasis->getMassMatrix(tempMass);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMass, massMatrix);
    massMatrixInv = massMatrix.inverse();

    // Store grad stiffness matrix in each direction
    std::vector<Eigen::MatrixXd> gradStiffnessMatrix(updateDirs.size());
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];
      Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, tempMatrix);
      gradStiffnessMatrix[dir] = Eigen::MatrixXd(nlocal, nlocal);

      copyLuceeToEigen(tempMatrix, gradStiffnessMatrix[dir]);
    }

    // get number of surface quadrature points
    int nSurfQuad = nodalBasis->getNumSurfGaussNodes();
    // get data needed for Gaussian quadrature
    int nVolQuad = nodalBasis->getNumGaussNodes();
    std::vector<double> volWeights(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NC);
    volQuad.reset(nVolQuad, nlocal, NC);

    double computationalScale = 1.0;
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      computationalScale *= grid.getDx(dimIndex);

    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex]/computationalScale;
    
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
    copyLuceeToEigen(tempVolCoords, volQuad.coordMat);

    std::vector<Eigen::MatrixXd> derivMatrices(updateDirs.size());

    // Compute gradients of basis functions evaluated at volume quadrature points
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];
      // Each row is a quadrature point; each column is a basis function with derivative applied
      // grid.getDx(dir) multiplied on to remove grid scale. Proper value needs to be multiplied in update()
      Eigen::MatrixXd derivMatrix = grid.getDx(dir)*volQuad.interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

      // Store gradients of the same basis function at various quadrature points
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        volQuad.pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);

      derivMatrices[d] = derivMatrix;
    }

    // Get data for surface quadrature
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];
      // temporary variables
      std::vector<double> tempSurfWeights(nSurfQuad);
      Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal);
      Lucee::Matrix<double> tempSurfCoords(nSurfQuad, NC);

      // Reset surface quadrature structures
      surfLowerQuad[dir].reset(nSurfQuad, nlocal, NC);
      surfUpperQuad[dir].reset(nSurfQuad, nlocal, NC);
      
      // Compute scale factor
      double surfAreaComputational = 1.0;
      for (int dirIndex = 0; dirIndex < NDIM; dirIndex++)
      {
        if (dirIndex != dir)
          surfAreaComputational = surfAreaComputational*grid.getDx(dirIndex);
      }

      // lower surface data
      nodalBasis->getSurfLowerGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfLowerQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex]/surfAreaComputational;
      copyLuceeToEigen(tempSurfQuad, surfLowerQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfLowerQuad[dir].coordMat);

      // upper surface data
      nodalBasis->getSurfUpperGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfUpperQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex]/surfAreaComputational;
      copyLuceeToEigen(tempSurfQuad, surfUpperQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfUpperQuad[dir].coordMat);
    }

    // Compute gradients of basis functions evaluated at surface quadrature points
    // surf keeps track of surfaces to evaluate derivatives on
    for (int s = 0; s < updateDirs.size(); s++)
    {
      int surf = updateDirs[s];
      for (int d = 0; d < updateDirs.size(); d++)
      {
        // dir is direction gradient is taken in
        int dir = updateDirs[d];
        // Each row is a quadrature point; each column is a basis function with derivative applied
        // grid.getDx(dir) multiplied on to remove grid scale. Proper value needs to be multiplied in update()
        Eigen::MatrixXd derivMatrix = grid.getDx(dir)*surfUpperQuad[surf].interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

        // Use derivMatrix to create matrices for gradients of basis functions
        for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
          surfUpperQuad[surf].pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);

        // Do the same for the lower surface
        derivMatrix = grid.getDx(dir)*surfLowerQuad[surf].interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();
        for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
          surfLowerQuad[surf].pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);
      }
    }

    // CONSIDER: instead of allocating NDIM size vectors, only store those needed in updateDirs
    bigStoredUpperSurfMatrices.resize(updateDirs.size());
    bigStoredLowerSurfMatrices.resize(updateDirs.size());
    bigStoredVolMatrices.resize(updateDirs.size());

    // Store three matrices at each cell
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];

      bigStoredUpperSurfMatrices[d] = computationalScale*
        (1.0/jacobianFactor)*massMatrixInv*surfUpperQuad[dir].interpMat.transpose();
      bigStoredLowerSurfMatrices[d] = computationalScale*
        (1.0/jacobianFactor)*massMatrixInv*surfLowerQuad[dir].interpMat.transpose();
      bigStoredVolMatrices[d] = computationalScale*
        (1.0/jacobianFactor)*massMatrixInv*derivMatrices[d].transpose();
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  PoissonBracketSimpleUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& aCurr = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& hamil = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& aNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    double dt = t-this->getCurrTime();
    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion();

    double cfla = 0.0; // maximum CFL number used

    int nlocal = nodalBasis->getNumNodes();
    int nVolQuad = nodalBasis->getNumGaussNodes();
    int nSurfQuad = nodalBasis->getNumSurfGaussNodes();

    Lucee::ConstFieldPtr<double> aCurrPtr = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_l = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> aCurrPtr_r = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilPtr = hamil.createConstPtr();
    Lucee::FieldPtr<double> aNewPtr = aNew.createPtr();
    Lucee::FieldPtr<double> aNewPtr_l = aNew.createPtr();
    Lucee::FieldPtr<double> aNewPtr_r = aNew.createPtr();

    aNew = 0.0; // use aNew to store increment initially

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    // Used in volume integral loop
    Eigen::VectorXd fVec(nlocal);
    Eigen::MatrixXd alpha(NDIM, nVolQuad);
    Eigen::VectorXd cflaVec(NDIM);
    Eigen::VectorXd dtdqVec(NDIM);

    // Contributions from volume integrals
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      aCurr.setPtr(aCurrPtr, idx);
      hamil.setPtr(hamilPtr, idx);
      grid.setIndex(idx);

      // Compute gradient of hamiltonian
      Eigen::MatrixXd hamilDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nVolQuad);
      for (int i = 0; i < nlocal; i++)
        hamilDerivAtQuad.noalias() += hamilPtr[i]*volQuad.pDiffMatrix[i];
      // Scale components based on grid
      for (int dirIndex = 0; dirIndex < NDIM; dirIndex++)
        hamilDerivAtQuad.row(dirIndex) *= grid.getSurfArea(dirIndex)/grid.getVolume();

      // Get alpha from appropriate function
      equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, volQuad.interpMat, idx, alpha);
      
      // Computing maximum cfla
      for (int i = 0; i < NDIM; i++)
        dtdqVec(i) = dt/(grid.getVolume()/grid.getSurfArea(i))/jacobianFactor;

      // Loop over each alpha vector evaluated at a quadrature point
      for (int i = 0; i < nVolQuad; i++)
      {
        double maxCflaInCol = alpha.col(i).cwiseAbs().cwiseProduct(dtdqVec).maxCoeff();
        if (maxCflaInCol > cfla)
          cfla = maxCflaInCol;
      }

      // Get a vector of f at quad points
      for (int i = 0; i < nlocal; i++)
        fVec(i) = aCurrPtr[i];
      Eigen::VectorXd fAtQuad = volQuad.interpMat*fVec;

      aNew.setPtr(aNewPtr, idx);
      for (int d = 0;  d < updateDirs.size(); d++)
      {
        int dir = updateDirs[d];
        // Scale factor
        double gridScale = grid.getSurfArea(dir)/grid.getVolume();
        for (int i = 0; i < nlocal; i++)
        {
          for (int qp = 0; qp < nVolQuad; qp++)
            aNewPtr[i] += gridScale*volQuad.weights(qp)*bigStoredVolMatrices[d](i,qp)*alpha(dir,qp)*fAtQuad(qp);
        }
      }
    }
     
    // Check to see if we need to retake time step
    if (cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);
     
    // Used in surface integral loop
    Eigen::VectorXd rightData(nlocal);
    Eigen::VectorXd leftData(nlocal);
    Eigen::MatrixXd alphaSurf(NDIM, nSurfQuad);
    Eigen::VectorXd numericalFluxAtQuad(nSurfQuad);

    // Contributions from surface integrals
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];
      // create sequencer to loop over *each* NDIM-1 slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seqLowerDim(localRgn.deflate(dir));
      // lower and upper bounds of NDIM-1 slice. (We need to make sure that flux
      // is computed for one edge outside domain interior)
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

      // Check to see if we have zero flux BCs in this direction
      // If true, then increase/decrease edge index by 1
      if (zeroFluxFlags[dir] == true)
      {
        if (sliceLower == globalRgn.getLower(dir))
          sliceLower = globalRgn.getLower(dir)+1;
        if (sliceUpper == globalRgn.getUpper(dir)+1)
          sliceUpper = globalRgn.getUpper(dir);
      }

      int idxr[NDIM];
      int idxl[NDIM];

      // loop over each 1D slice
      // TODO: use only lower/upper or right/left naming convention
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idxr);
        seqLowerDim.fillWithIndex(idxl);
        // loop over each edge
        for (int sliceIndex = sliceLower; sliceIndex < sliceUpper; sliceIndex++)
        { 
          idxr[dir] = sliceIndex;
          idxl[dir] = sliceIndex-1;

          aCurr.setPtr(aCurrPtr_r, idxr);
          aCurr.setPtr(aCurrPtr_l, idxl);
          // Copy data to Eigen vectors
          for (int i = 0; i < nlocal; i++)
          {
              rightData(i) = aCurrPtr_r[i];
              leftData(i) = aCurrPtr_l[i];
          }

          // Compute alpha at edge quadrature nodes (making use of alpha dot n being continuous)
          if (sliceIndex != globalRgn.getUpper(dir))
          {
            // Hamiltonian is continuous, so use right always
            // UNLESS left cell is globalRgn.getUpper(dir)
            hamil.setPtr(hamilPtr, idxr);
            grid.setIndex(idxr);
            // Compute gradient of hamiltonian at surface nodes
            Eigen::MatrixXd hamilDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nSurfQuad);
            for (int i = 0; i < nlocal; i++)
              hamilDerivAtQuad.noalias() += hamilPtr[i]*surfLowerQuad[dir].pDiffMatrix[i];
            // Scale components based on grid
            for (int dirIndex = 0; dirIndex < NDIM; dirIndex++)
              hamilDerivAtQuad.row(dirIndex) *= grid.getSurfArea(dirIndex)/grid.getVolume();
            equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfLowerQuad[dir].interpMat, idxr, alphaSurf);
          }
          else
          {
            // Hamiltonian is continuous, so use right always
            // UNLESS left cell is globalRgn.getUpper(dir)
            hamil.setPtr(hamilPtr, idxl);
            grid.setIndex(idxl);
            // Compute gradient of hamiltonian at surface nodes
            Eigen::MatrixXd hamilDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nSurfQuad);
            for (int i = 0; i < nlocal; i++)
              hamilDerivAtQuad.noalias() += hamilPtr[i]*surfUpperQuad[dir].pDiffMatrix[i];
            // Scale components based on grid
            for (int dirIndex = 0; dirIndex < NDIM; dirIndex++)
              hamilDerivAtQuad.row(dirIndex) *= grid.getSurfArea(dirIndex)/grid.getVolume();
            equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfUpperQuad[dir].interpMat, idxl, alphaSurf);
          }

          // Construct normal vector
          Eigen::VectorXd normalVec = Eigen::VectorXd::Zero(NDIM);
          normalVec(dir) = 1.0;
          
          // Calculate alphaDotN at all quadrature points
          Eigen::RowVectorXd alphaDotN = normalVec.transpose()*alphaSurf;

          // Compute numerical flux
          computeNumericalFlux(alphaDotN, surfUpperQuad[dir].interpMat*leftData,
              surfLowerQuad[dir].interpMat*rightData, numericalFluxAtQuad);

          aNew.setPtr(aNewPtr_r, idxr);
          aNew.setPtr(aNewPtr_l, idxl);
          // Surface area will be the same for both cells
          grid.setIndex(idxl);
          double gridScaleLeft = grid.getSurfArea(dir)/grid.getVolume();
          grid.setIndex(idxr);
          double gridScaleRight = grid.getSurfArea(dir)/grid.getVolume();
          for (int i = 0; i < nlocal; i++)
          {
            for (int qp = 0; qp < nSurfQuad; qp++)
            {
              aNewPtr_l[i] -= bigStoredUpperSurfMatrices[d](i,qp)*
                numericalFluxAtQuad(qp)*surfUpperQuad[dir].weights(qp)*gridScaleLeft;
              aNewPtr_r[i] += bigStoredLowerSurfMatrices[d](i,qp)*
                numericalFluxAtQuad(qp)*surfLowerQuad[dir].weights(qp)*gridScaleRight;
            }
          }
        }
      }
    }
    

    // NOTE: If only calculation of increments are requested, the final
    // update is not performed. This means that the multiplication
    // of the DG RHS with dt is not done, something to keep in mind if
    // using the increment in time-dependent update.
    seq = Lucee::RowMajorSequencer<NDIM>(localRgn);

    if (onlyIncrement == false)
    {
      while(seq.step())
      {
        seq.fillWithIndex(idx);
        aCurr.setPtr(aCurrPtr, idx);
        aNew.setPtr(aNewPtr, idx);

        for (int i = 0; i < nlocal; i++)
          aNewPtr[i] = aCurrPtr[i] + dt*aNewPtr[i];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }
  
  template <unsigned NDIM>
  void
  PoissonBracketSimpleUpdater<NDIM>::declareTypes()
  {
    // takes two inputs (aCurr, hamil)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns one output, (aNew)
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  PoissonBracketSimpleUpdater<NDIM>::computeNumericalFlux(const Eigen::VectorXd& alphaDotN,
    const Eigen::VectorXd& leftValsAtQuad, const Eigen::VectorXd& rightValsAtQuad,
    Eigen::VectorXd& numericalFluxAtQuad)
  {
    double leftVal;
    double rightVal;
    // Loop through all quadrature points
    for (int quadIndex = 0; quadIndex < numericalFluxAtQuad.size(); quadIndex++)
    {
      //set distribution function at quad to 0 if value is negative
      //if (leftValsAtQuad(quadIndex) < 0.0)
      //  leftVal = 0.0;
      //else
        leftVal = leftValsAtQuad(quadIndex);

      //if (rightValsAtQuad(quadIndex) < 0.0)
      //  rightVal = 0.0;
      //else
        rightVal = rightValsAtQuad(quadIndex);
      
      if (fluxType == UPWIND)
      {
        if (alphaDotN(quadIndex) > 0.0)
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*leftVal;
        else
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*rightVal;
      }
      else if (fluxType == CENTRAL)
        numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*0.5*
          (rightVal + leftVal);
      else if (fluxType == DOWNWIND)
      {
        if (alphaDotN(quadIndex) > 0.0)
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*rightVal;
        else
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*leftVal;
      }
    }
  }

  template <unsigned NDIM>
  void
  PoissonBracketSimpleUpdater<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class PoissonBracketSimpleUpdater<1>;
  template class PoissonBracketSimpleUpdater<2>;
  template class PoissonBracketSimpleUpdater<3>;
  template class PoissonBracketSimpleUpdater<4>;
  template class PoissonBracketSimpleUpdater<5>;
}

/**
 * @file	LcPoissonBracketImpUpdater.cpp
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 * This version is used to test speed improvements to the existing updater.
 * NOTE: Jacobian must be supplied
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMathLib.h>
#include <LcPoissonBracketImpUpdater.h>

#include <ctime>

namespace Lucee
{
  static const unsigned UPWIND = 0;
  static const unsigned CENTRAL = 1;
  static const unsigned DOWNWIND = 2;

// set id for module system
  template <> const char *PoissonBracketImpUpdater<1>::id = "PoissonBracketImp1D";
  template <> const char *PoissonBracketImpUpdater<2>::id = "PoissonBracketImp2D";
  template <> const char *PoissonBracketImpUpdater<3>::id = "PoissonBracketImp3D";
  template <> const char *PoissonBracketImpUpdater<4>::id = "PoissonBracketImp4D";
  template <> const char *PoissonBracketImpUpdater<5>::id = "PoissonBracketImp5D";

  template <unsigned NDIM>
  PoissonBracketImpUpdater<NDIM>::PoissonBracketImpUpdater()
    : UpdaterIfc()
  {
    totalVolTime = totalSurfTime = jacAtQuad = 0.0;
    computeVolAlphaAtQuadNodesTime = computeSurfAlphaAtQuadNodesTime = 0.0;
    vol_loop1 = vol_loop2 = vol_loop3 = vol_loop4 = 0.0;
    surf_loop1 = surf_loop2 = surf_loop3 = surf_loop4 = 0.0;
  }

  template <unsigned NDIM>
  PoissonBracketImpUpdater<NDIM>::~PoissonBracketImpUpdater()
  {/*
    std::cout << "Total volume integral time = " << totalVolTime << std::endl;
    std::cout << "Total surface integral time = " << totalSurfTime << std::endl;
    std::cout << "Total alpha for volume integral time = " << computeVolAlphaAtQuadNodesTime << std::endl;
    std::cout << "Total alpha for surface integral time = " << computeSurfAlphaAtQuadNodesTime << std::endl;
    std::cout << "Vol loops  = "
              << vol_loop1 << " + "
              << vol_loop2 << " + "
              << vol_loop3 << " + "
              << vol_loop4 << " = "
              << (vol_loop1+vol_loop2+vol_loop3+vol_loop4)
              << std::endl;

    std::cout << "Surf loops  = "
              << surf_loop1 << " + "
              << surf_loop2 << " + "
              << surf_loop3 << " + "
              << surf_loop4 << " = "
              << (surf_loop1+surf_loop2+surf_loop3+surf_loop4)
              << std::endl;*/
  }  
  
  template <unsigned NDIM>
  void 
  PoissonBracketImpUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("PoissonBracketImpUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasObject<Lucee::PoissonBracketEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::PoissonBracketEquation>("equation");
    else
    {
      Lucee::Except lce("PoissonBracketImpUpdater::readInput: Must specify an equation to solve!");
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
        Lucee::Except lce("PoissonBracketImpUpdater::readInput: 'fluxType' ");
        lce << tbl.getString("fluxType") << " is not valid";
        throw lce;
      }
    }

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    hasJacobian = false;
    if (tbl.hasObject<Lucee::Field<NDIM, double> >("jacobianField"))
    {
      hasJacobian = true;
      jacobianField = &tbl.getObject<Lucee::Field<NDIM, double> >("jacobianField");
    }
    else
    {
      throw Lucee::Except("PoissonBracketImpUpdater::readInput: A jacobianField MUST be supplied!");
    }

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
          throw Lucee::Except("PoissonBracketImpUpdater::readInput: updateDirections must be a table with each element < NDIM");
      }
    }
    else
    {
      for (int i = 0; i < NDIM; i++)
        updateDirs.push_back(i);
    }

    zeroFluxFlagsLower = std::vector<bool>(NDIM);
    zeroFluxFlagsUpper = std::vector<bool>(NDIM);
    // Default is no zero flux bcs in any directions
    for (int i = 0; i < NDIM; i++)
    {
      zeroFluxFlagsLower[i] = false;
      zeroFluxFlagsUpper[i] = false;
    }

    if (tbl.hasNumVec("zeroFluxDirectionsUpper"))
    {
      std::vector<double> tempDirs = tbl.getNumVec("zeroFluxDirectionsUpper");
      for (int i = 0; i < std::min<int>(NDIM, tempDirs.size()); i++)
      {
        int d = (int) tempDirs[i];

        if (d < NDIM)
          zeroFluxFlagsUpper[d] = true;
        else
          throw Lucee::Except("PoissonBracketImpUpdater::readInput: zeroFluxDirectionsUpper must be a table with each element < NDIM");
      }
    }

    if (tbl.hasNumVec("zeroFluxDirectionsLower"))
    {
      std::vector<double> tempDirs = tbl.getNumVec("zeroFluxDirectionsLower");
      for (int i = 0; i < std::min<int>(NDIM, tempDirs.size()); i++)
      {
        int d = (int) tempDirs[i];

        if (d < NDIM)
          zeroFluxFlagsLower[d] = true;
        else
          throw Lucee::Except("PoissonBracketImpUpdater::readInput: zeroFluxDirectionsLower must be a table with each element < NDIM");
      }
    }
  }

  template <unsigned NDIM>
  void 
  PoissonBracketImpUpdater<NDIM>::initialize()
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
      gradStiffnessMatrix[d] = Eigen::MatrixXd(nlocal, nlocal);

      copyLuceeToEigen(tempMatrix, gradStiffnessMatrix[d]);
    }

    // get number of surface quadrature points
    int nSurfQuad = nodalBasis->getNumSurfGaussNodes();
    // get data needed for Gaussian quadrature
    int nVolQuad = nodalBasis->getNumGaussNodes();
    std::vector<double> volWeights(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NC);
    volQuad.reset(nVolQuad, nlocal, NC);

    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex];
    
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
    copyLuceeToEigen(tempVolCoords, volQuad.coordMat);

    std::vector<Eigen::MatrixXd> derivMatrices(updateDirs.size());

    // Get data for surface quadrature
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];
      // temporary variables
      std::vector<double> tempSurfWeights(nSurfQuad);
      Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal);
      Lucee::Matrix<double> tempSurfCoords(nSurfQuad, NC);

      int nSurfNodes = nodalBasis->getNumSurfLowerNodes(dir);

      // Reset surface quadrature structures
      surfLowerQuad[dir].reset(nSurfQuad, nSurfNodes, NC);
      surfUpperQuad[dir].reset(nSurfQuad, nSurfNodes, NC);
      
      // lower surface data
      nodalBasis->getSurfLowerGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      nodalBasis->getSurfLowerNodeNums(dir, surfLowerQuad[dir].nodeNums);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
      {
        surfLowerQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex];
        for (int nodeIndex = 0; nodeIndex < nSurfNodes; nodeIndex++)
          surfLowerQuad[dir].interpMat(quadIndex, nodeIndex) = 
            tempSurfQuad(quadIndex, surfLowerQuad[dir].nodeNums[nodeIndex]);
      }
      //copyLuceeToEigen(tempSurfQuad, surfLowerQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfLowerQuad[dir].coordMat);

      // upper surface data
      nodalBasis->getSurfUpperGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      nodalBasis->getSurfUpperNodeNums(dir, surfUpperQuad[dir].nodeNums);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
      {
        surfUpperQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex];
        for (int nodeIndex = 0; nodeIndex < nSurfNodes; nodeIndex++)
          surfUpperQuad[dir].interpMat(quadIndex, nodeIndex) = 
            tempSurfQuad(quadIndex, surfUpperQuad[dir].nodeNums[nodeIndex]);
      }
      //copyLuceeToEigen(tempSurfQuad, surfUpperQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfUpperQuad[dir].coordMat);
    }

    // IMPROVED LOOP
    gradMatrices.resize(updateDirs.size());
    for (int d = 0; d < updateDirs.size(); d++)
    {
      // column 'k' of gradMatrices[d] contains the basis set expansion
      // of the dth direction derivative of basis function k
      gradMatrices[d] = massMatrixInv*gradStiffnessMatrix[d].transpose();
      // Each row is a quadrature point, each column is a basis function derivative
      derivMatrices[d] = volQuad.interpMat*gradMatrices[d];
    }

    // TotalSize keeps track of total number of cells we will need to store matrices for
    // (this is physical domain of region + 1 layer of ghost cells)
    int totalSize = 1;
    for (int dir = 0; dir < NDIM; dir++)
    {
      localRgn.setLower(dir, localRgn.getLower(dir)-1);
      localRgn.setUpper(dir, localRgn.getUpper(dir)+1);
      totalSize *= localRgn.getUpper(dir) - localRgn.getLower(dir);
    }
    
    Lucee::RowMajorSequencer<NDIM> volSeq = RowMajorSequencer<NDIM>(localRgn);
    Lucee::RowMajorIndexer<NDIM> volIdxr = RowMajorIndexer<NDIM>(localRgn);
    
    // CONSIDER: instead of allocating NDIM size vectors, only store those needed in updateDirs
    bigStoredUpperSurfMatrices.resize(totalSize, std::vector<Eigen::MatrixXd>(updateDirs.size()));
    bigStoredLowerSurfMatrices.resize(totalSize, std::vector<Eigen::MatrixXd>(updateDirs.size()));
    bigStoredVolMatrices.resize(totalSize, std::vector<Eigen::MatrixXd>(updateDirs.size()));

    Lucee::ConstFieldPtr<double> jacobianPtr = jacobianField->createConstPtr();

    while(volSeq.step())
    {
      volSeq.fillWithIndex(idx);
      int cellIndex = volIdxr.getIndex(idx);

      // Compute Jacobian-weighted mass matrix inverse
      if (hasJacobian == true)
      {
        jacobianField->setPtr(jacobianPtr, idx);
        Eigen::VectorXd jacobianVec(nlocal);
        for (int i = 0; i < nlocal; i++)
          jacobianVec(i) = jacobianPtr[i];

        Eigen::MatrixXd tempJMassMatrix(nlocal, nlocal);

        Eigen::VectorXd jacobianAtQuad = volQuad.interpMat*jacobianVec;

        // Use jacobian to compute new massMatrixInv
        for (int j = 0; j < nlocal; j++)
          for (int k = 0; k < nlocal; k++)
            tempJMassMatrix(j,k) = volQuad.interpMat.col(j).cwiseProduct(volQuad.interpMat.col(k)).cwiseProduct(jacobianAtQuad).dot(volQuad.weights);

        massMatrixInv = tempJMassMatrix.inverse();
      }

      // Store three matrices at each cell
      for (int d = 0; d < updateDirs.size(); d++)
      {
        int dir = updateDirs[d];

        // temporary variables
        std::vector<double> tempSurfWeights(nSurfQuad);
        Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal);
        Lucee::Matrix<double> tempSurfCoords(nSurfQuad, NC);

        Eigen::MatrixXd tempInterpMatLower(nSurfQuad, nlocal);
        nodalBasis->getSurfLowerGaussQuadData(dir, tempSurfQuad,
          tempSurfCoords, tempSurfWeights);
        copyLuceeToEigen(tempSurfQuad, tempInterpMatLower);
        
        Eigen::MatrixXd tempInterpMatUpper(nSurfQuad, nlocal);
        nodalBasis->getSurfUpperGaussQuadData(dir, tempSurfQuad,
          tempSurfCoords, tempSurfWeights);
        copyLuceeToEigen(tempSurfQuad, tempInterpMatUpper);

        bigStoredUpperSurfMatrices[cellIndex][d] = massMatrixInv*tempInterpMatUpper.transpose();
        bigStoredLowerSurfMatrices[cellIndex][d] = massMatrixInv*tempInterpMatLower.transpose();
        // Multiply each row of the surface matrices by quadrature weights so we only need to multiply
        // this matrix with the numerical flux vector to get update weights
        for (int rowIndex = 0; rowIndex < bigStoredUpperSurfMatrices[cellIndex][d].rows(); rowIndex++)
        {
          bigStoredUpperSurfMatrices[cellIndex][d].row(rowIndex) = 
            bigStoredUpperSurfMatrices[cellIndex][d].row(rowIndex).
            cwiseProduct(surfUpperQuad[dir].weights.transpose());
          bigStoredLowerSurfMatrices[cellIndex][d].row(rowIndex) = 
            bigStoredLowerSurfMatrices[cellIndex][d].row(rowIndex).
            cwiseProduct(surfLowerQuad[dir].weights.transpose());
        }
        // Each row is a basis function derivative, each column is a quadrature point location
        bigStoredVolMatrices[cellIndex][d] = massMatrixInv*derivMatrices[d].transpose();
      }
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  PoissonBracketImpUpdater<NDIM>::update(double t)
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

    Lucee::Region<NDIM, int> effectiveRgn(localRgn);
    for (int dir = 0; dir < NDIM; dir++)
    {
      effectiveRgn.setLower(dir, localRgn.getLower(dir)-1);
      effectiveRgn.setUpper(dir, localRgn.getUpper(dir)+1);
    }

    Lucee::RowMajorIndexer<NDIM> volIdxr = RowMajorIndexer<NDIM>(effectiveRgn);

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
    // Testing: see LcWavePropagationUpdater
    Lucee::ConstFieldPtr<double> jacobianPtr = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> jacobianPtr_l = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> jacobianPtr_r = aCurr.createConstPtr();

    aNew = 0.0; // use aNew to store increment initially

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    // Potentially costly Eigen structures for surface integral calculations  
    // Components of solution on left and right of boundary
    Eigen::MatrixXd hamilDerivAtQuad = Eigen::MatrixXd(NDIM, nSurfQuad);
    Eigen::MatrixXd hamilDeriv = Eigen::MatrixXd::Zero(NDIM, nlocal);
    Eigen::MatrixXd alpha(NDIM, nVolQuad);
    Eigen::RowVectorXd alphaDotN(nSurfQuad);
    Eigen::MatrixXd alphaSurf(NDIM, nSurfQuad);
    Eigen::VectorXd normalVec = Eigen::VectorXd(NDIM);
    Eigen::VectorXd numericalFluxAtQuad(nSurfQuad);
    Eigen::VectorXd upperResultVector(nlocal);
    Eigen::VectorXd lowerResultVector(nlocal);
    Eigen::VectorXd jacobianVec(nlocal);
    Eigen::VectorXd hamilVec(nlocal);
    Eigen::VectorXd jacobianAtQuad(nVolQuad);
    Eigen::MatrixXd hamilDerivAtQuadVol(NDIM, nVolQuad);
    Eigen::VectorXd fVec(nlocal);
    Eigen::VectorXd resultVector(nlocal);
    Eigen::VectorXd fAtQuad(nVolQuad);

    bool needToRetakeStep = false;

    // Contributions from volume integrals
    while(seq.step())
    {
      clock_t tm_totalVolTime_s = clock();
      seq.fillWithIndex(idx);
      aCurr.setPtr(aCurrPtr, idx);
      hamil.setPtr(hamilPtr, idx);
      int cellIndex = volIdxr.getIndex(idx);

      jacobianAtQuad = Eigen::VectorXd::Ones(nVolQuad);
      // Figure out Jacobian at quadrature points so its effect can be
      // accounted for when checking CFL condition
      if (hasJacobian == true)
      {
        jacobianField->setPtr(jacobianPtr, idx);
        for (int i = 0; i < nlocal; i++)
          jacobianVec(i) = jacobianPtr[i];

        jacobianAtQuad = volQuad.interpMat*jacobianVec;
      }

      // Store hamiltonian as eigen vector
      for (int i = 0; i < nlocal; i++)
        hamilVec(i) = hamilPtr[i];
      // Loop through all directions
      // Compute gradient of hamiltonian
      clock_t tm_1_s = clock();
      for (int dir = 0; dir < updateDirs.size(); dir++)
        hamilDeriv.row(updateDirs[dir]) = gradMatrices[dir]*hamilVec;
      vol_loop1 += (double) (clock() - tm_1_s)/CLOCKS_PER_SEC;
      // Compute gradient of hamiltonian at all volume quadrature points
      clock_t tm_2_s = clock();
      hamilDerivAtQuadVol.noalias() = hamilDeriv*volQuad.interpMat.transpose();
      vol_loop2 += (double) (clock() - tm_2_s)/CLOCKS_PER_SEC;

      clock_t tm_computeAlphaAtQuadNodesTime_s = clock();
      // Get alpha from appropriate function
      equation->computeAlphaAtQuadNodes(hamilDerivAtQuadVol, volQuad.interpMat, idx, alpha);
      clock_t tm_computeAlphaAtQuadNodesTime_e = clock();
      computeVolAlphaAtQuadNodesTime +=
        (double) (tm_computeAlphaAtQuadNodesTime_e-tm_computeAlphaAtQuadNodesTime_s)/CLOCKS_PER_SEC;
      
      // Computing maximum cfla
      grid.setIndex(idx);
      Eigen::VectorXd dtdqVec(NDIM);
      for (int i = 0; i < NDIM; i++)
        dtdqVec(i) = dt/grid.getDx(i);

      clock_t tm_3_s = clock();
      Eigen::VectorXd cflaVec(NDIM);
      // Loop over each alpha vector evaluated at a quadrature point
      for (int i = 0; i < nVolQuad; i++)
      {
        double maxCflaInCol = alpha.col(i).cwiseAbs().cwiseProduct(dtdqVec).maxCoeff()/fabs(jacobianAtQuad(i));
        if (maxCflaInCol > cfla)
        {
          cfla = maxCflaInCol;
          if (cfla > cflm)
            needToRetakeStep = true;
        }
      }
      vol_loop3 += (double) (clock() - tm_3_s)/CLOCKS_PER_SEC;

      // Get a vector of f at quad points
      for (int i = 0; i < nlocal; i++)
        fVec(i) = aCurrPtr[i];
      // Compute distribution function at quadrature points cwise multiplied by quadrature weights
      fAtQuad = volQuad.weights.cwiseProduct(volQuad.interpMat*fVec);

      aNew.setPtr(aNewPtr, idx);
      clock_t tm_4_s = clock();
      for (int d = 0;  d < updateDirs.size(); d++)
      {
        int dir = updateDirs[d];
        resultVector = bigStoredVolMatrices[cellIndex][d]*alpha.row(dir).transpose().cwiseProduct(fAtQuad);

        for (int i = 0; i < nlocal; i++)
          aNewPtr[i] += resultVector(i);
      }
      vol_loop4 += (double) (clock() - tm_4_s)/CLOCKS_PER_SEC;

      clock_t tm_totalVolTime_e = clock();
      totalVolTime += (double) (tm_totalVolTime_e - tm_totalVolTime_s)/CLOCKS_PER_SEC;
      // Loop over surfaces on which integrals need to be calculated
      // Only do this if cfl limit has not been violated yet
      clock_t tm_totalSurfTime_s = clock();
      if (needToRetakeStep == false)
      {
        for (int d = 0; d < updateDirs.size(); d++)
        {
          int dir = updateDirs[d];
          int idxl[NDIM];
          int idxr[NDIM];
          
          for (int i = 0; i < NDIM; i++)
            idxl[i] = idx[i];
          idxl[dir] = idx[dir] - 1;

          // Don't need to do anything if zeroFluxFlag in this direction and on lower boundary
          if (idx[dir] == globalRgn.getLower(dir) && zeroFluxFlagsLower[dir] == true )
            continue;
          // Otherwise, need to compute numerical flux. Get data from lower element
          aCurr.setPtr(aCurrPtr_l, idxl);

          int nSurfNodes = nodalBasis->getNumSurfLowerNodes(dir);
          // Copy data to Eigen vectors
          Eigen::VectorXd leftData(nSurfNodes);
          Eigen::VectorXd rightData(nSurfNodes);
          for (int i = 0; i < nSurfNodes; i++)
          {
            leftData(i) = aCurrPtr_l[ surfUpperQuad[dir].nodeNums[i] ];
            rightData(i) = aCurrPtr[ surfLowerQuad[dir].nodeNums[i] ];
          }
          // Get linear indices
          int cellIndexLeft = volIdxr.getIndex(idxl);
          // Compute reduced hamiltonian derivative structure
          Eigen::MatrixXd hamilDerivSmall(NDIM, nSurfNodes);
          for (int i = 0; i < nSurfNodes; i++)
            hamilDerivSmall.col(i) = hamilDeriv.col( surfLowerQuad[dir].nodeNums[i] );
          // Compute gradient of hamiltonian at all (lower) surface quadrature points
          clock_t tm_1_s = clock();
          hamilDerivAtQuad.noalias() = hamilDerivSmall*surfLowerQuad[dir].interpMat.transpose();
          surf_loop1 += (double) (clock()-tm_1_s)/CLOCKS_PER_SEC;
          // Compute alpha at edge quadrature nodes (making use of alpha dot n being continuous)
          clock_t tm_computeAlphaAtQuadNodesTime_s = clock();    
          equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfLowerQuad[dir].interpMat,
              surfLowerQuad[dir].nodeNums, idx, alphaDotN, dir);
          clock_t tm_computeAlphaAtQuadNodesTime_e = clock();
          computeSurfAlphaAtQuadNodesTime +=
            (double) (tm_computeAlphaAtQuadNodesTime_e-tm_computeAlphaAtQuadNodesTime_s)/CLOCKS_PER_SEC;

          //equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfLowerQuad[dir].interpMat, idx, alphaSurf);
          // Construct normal vector
          //normalVec.setZero(NDIM);
          //normalVec(dir) = 1.0;
          // Calculate alphaDotN at all quadrature points
          //Eigen::RowVectorXd alphaDotNTemp = normalVec.transpose()*alphaSurf;

          // Compare alphaDotN and alphaDotNTemp
          //std::cout << "alpha a " << dir << std::endl << alphaDotN << std::endl << "alpha b " << dir << std::endl << alphaDotNTemp << std::endl;
          // TODO: account for jacobian factor, since it has been multiplied into alpha already?
          // Don't need to do this for now, as long as jacobian factor is > 0

          // Compute numerical flux
          computeNumericalFlux(alphaDotN, surfUpperQuad[dir].interpMat*leftData,
              surfLowerQuad[dir].interpMat*rightData, numericalFluxAtQuad);
          
          clock_t tm_2_s = clock();
          upperResultVector = bigStoredUpperSurfMatrices[cellIndexLeft][d]*numericalFluxAtQuad;
          lowerResultVector = bigStoredLowerSurfMatrices[cellIndex][d]*numericalFluxAtQuad;
          surf_loop2 += (double) (clock()-tm_2_s)/CLOCKS_PER_SEC;

          clock_t tm_3_s = clock();
          // Accumulate to solution
          aNew.setPtr(aNewPtr_r, idx);
          aNew.setPtr(aNewPtr_l, idxl);
          for (int i = 0; i < nlocal; i++)
          {
            aNewPtr_l[i] -= upperResultVector(i);
            aNewPtr_r[i] += lowerResultVector(i);
          }
          surf_loop3 += (double) (clock()-tm_3_s)/CLOCKS_PER_SEC;

          clock_t tm_4_s = clock();
          // Need to do integrals on the upper surface if we are on the upper boundary
          // because volume sequencer does not loop over ghost cells
          if(idx[dir] == localRgn.getUpper(dir)-1)
          {
            // Don't need to do anything only if zeroFluxFlags in this direction and on global boundary
            if (idx[dir] == globalRgn.getUpper(dir)-1 && zeroFluxFlagsUpper[dir] == true)
              continue;

            for (int i = 0; i < NDIM; i++)
            {
              idxr[i] = idx[i];
              idxl[i] = idx[i];
            }
            // Now, the right element is in the ghost region
            idxr[dir] = idx[dir] + 1;
            aCurr.setPtr(aCurrPtr_l, idxl);
            aCurr.setPtr(aCurrPtr_r, idxr);
            // Copy data to Eigen vectors
            for (int i = 0; i < nSurfNodes; i++)
            {
              leftData(i) = aCurrPtr_l[ surfUpperQuad[dir].nodeNums[i] ];
              rightData(i) = aCurrPtr_r[ surfLowerQuad[dir].nodeNums[i] ];
            }
            int cellIndexRight = volIdxr.getIndex(idxr);
            cellIndexLeft = volIdxr.getIndex(idxl);
            // Compute reduced hamiltonian derivative structure
            for (int i = 0; i < nSurfNodes; i++)
              hamilDerivSmall.col(i) = hamilDeriv.col( surfUpperQuad[dir].nodeNums[i] );

            hamilDerivAtQuad.noalias() = hamilDerivSmall*surfUpperQuad[dir].interpMat.transpose();
            
            // Compute alpha at edge quadrature nodes (making use of alpha dot n being continuous)
            equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfUpperQuad[dir].interpMat, 
                surfUpperQuad[dir].nodeNums, idxl, alphaDotN, dir);
            
            // Calculate alphaDotN at all quadrature points
            //alphaDotN = normalVec.transpose()*alphaSurf;
            // TODO: account for jacobian factor, since it has been multiplied into alpha already?
            // Don't need to do this for now, as long as jacobian factor is > 0

            // Compute numerical flux
            computeNumericalFlux(alphaDotN, surfUpperQuad[dir].interpMat*leftData,
                surfLowerQuad[dir].interpMat*rightData, numericalFluxAtQuad);
            
            upperResultVector = bigStoredUpperSurfMatrices[cellIndexLeft][d]*
              numericalFluxAtQuad;
            lowerResultVector = bigStoredLowerSurfMatrices[cellIndexRight][d]*
              numericalFluxAtQuad;

            // Accumulate to solution
            aNew.setPtr(aNewPtr_r, idxr);
            aNew.setPtr(aNewPtr_l, idxl);
            for (int i = 0; i < nlocal; i++)
            {
              aNewPtr_l[i] -= upperResultVector(i);
              aNewPtr_r[i] += lowerResultVector(i);
            }
          }
          surf_loop4 += (double) (clock()-tm_4_s)/CLOCKS_PER_SEC;
        }
      }
      clock_t tm_totalSurfTime_e = clock();
      totalSurfTime += (double) (tm_totalSurfTime_e - tm_totalSurfTime_s)/CLOCKS_PER_SEC;
    }
    
    // Check to see if we need to retake time step
    if (needToRetakeStep == true)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);
     
    
    /*
    clock_t tm_totalSurfTime_s = clock();
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

      // Loop over each 1D slice
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
          
          int cellIndexRight = volIdxr.getIndex(idxr);
          int cellIndexLeft = volIdxr.getIndex(idxl);
          // Hamiltonian is continuous, so use right value always
          hamil.setPtr(hamilPtr, idxr);
          // Copy data to Eigen vectors
          for (int i = 0; i < nlocal; i++)
          {
            rightData(i) = aCurrPtr_r[i];
            leftData(i) = aCurrPtr_l[i];
          }

          clock_t tm_1_s = clock();
          // Compute gradient of hamiltonian in right element
          //Eigen::MatrixXd hamilDeriv = Eigen::MatrixXd(NDIM, nlocal);
          // Store hamiltonian as eigen vector
          Eigen::VectorXd hamilVec(nlocal);
          for (int i = 0; i < nlocal; i++)
            hamilVec(i) = hamilPtr[i];
          // Loop through all directions
          for (int j = 0; j < updateDirs.size(); j++)
            hamilDeriv.row(updateDirs[j]) = gradMatrices[j]*hamilVec;
          // Compute gradient of hamiltonian at all surface quadrature points
          Eigen::MatrixXd hamilDerivAtQuad = hamilDeriv*surfLowerQuad[dir].interpMat.transpose();
          
          surf_loop1 += (double) (clock()-tm_1_s)/CLOCKS_PER_SEC;

          clock_t tm_computeAlphaAtQuadNodesTime_s = clock();          
          // Compute alpha at edge quadrature nodes (making use of alpha dot n being continuous)
          equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfLowerQuad[dir].interpMat, idxr, alpha);
          clock_t tm_computeAlphaAtQuadNodesTime_e = clock();
          computeSurfAlphaAtQuadNodesTime +=
            (double) (tm_computeAlphaAtQuadNodesTime_e-tm_computeAlphaAtQuadNodesTime_s)/CLOCKS_PER_SEC;

          // Construct normal vector
          normalVec.setZero(NDIM);
          normalVec(dir) = 1.0;
          
          // Calculate alphaDotN at all quadrature points
          Eigen::RowVectorXd alphaDotN = normalVec.transpose()*alpha;
          // TODO: account for jacobian factor, since it has been multiplied into alpha already?
          // Don't need to do this for now, as long as jacobian factor is > 0

          clock_t tm_2_s = clock();
          // Compute numerical flux
          computeNumericalFlux(alphaDotN, surfUpperQuad[dir].interpMat*leftData,
              surfLowerQuad[dir].interpMat*rightData, numericalFluxAtQuad);
          surf_loop2 += (double) (clock()-tm_2_s)/CLOCKS_PER_SEC;
          
          clock_t tm_3_s = clock();          
          Eigen::VectorXd upperResultVector = bigStoredUpperSurfMatrices[cellIndexLeft][d]*
            numericalFluxAtQuad.cwiseProduct(surfUpperQuad[dir].weights);
          Eigen::VectorXd lowerResultVector = bigStoredLowerSurfMatrices[cellIndexRight][d]*
            numericalFluxAtQuad.cwiseProduct(surfLowerQuad[dir].weights);
          surf_loop3 += (double) (clock()-tm_3_s)/CLOCKS_PER_SEC;

          // Accumulate to solution
          aNew.setPtr(aNewPtr_r, idxr);
          aNew.setPtr(aNewPtr_l, idxl);
          for (int i = 0; i < nlocal; i++)
          {
            aNewPtr_l[i] -= upperResultVector(i);
            aNewPtr_r[i] += lowerResultVector(i);
          }
        }
      }
    }
    clock_t tm_totalSurfTime_e = clock();
    totalSurfTime += (double) (tm_totalSurfTime_e - tm_totalSurfTime_s)/CLOCKS_PER_SEC;
    */

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
  PoissonBracketImpUpdater<NDIM>::declareTypes()
  {
    // takes two inputs (aCurr, hamil)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns one output, (aNew)
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  PoissonBracketImpUpdater<NDIM>::computeNumericalFlux(const Eigen::VectorXd& alphaDotN,
    const Eigen::VectorXd& leftValsAtQuad, const Eigen::VectorXd& rightValsAtQuad,
    Eigen::VectorXd& numericalFluxAtQuad)
  {
    // Loop through all quadrature points
    for (int quadIndex = 0; quadIndex < numericalFluxAtQuad.size(); quadIndex++)
    {
      if (fluxType == UPWIND)
      {
        if (alphaDotN(quadIndex) > 0.0)
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*leftValsAtQuad(quadIndex);
        else
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*rightValsAtQuad(quadIndex);
      }
      else if (fluxType == CENTRAL)
        numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*0.5*
          (rightValsAtQuad(quadIndex) + leftValsAtQuad(quadIndex));
      else if (fluxType == DOWNWIND)
      {
        if (alphaDotN(quadIndex) > 0.0)
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*rightValsAtQuad(quadIndex);
        else
          numericalFluxAtQuad(quadIndex) = alphaDotN(quadIndex)*leftValsAtQuad(quadIndex);
      }
    }
  }

  template <unsigned NDIM>
  void
  PoissonBracketImpUpdater<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class PoissonBracketImpUpdater<1>;
  template class PoissonBracketImpUpdater<2>;
  template class PoissonBracketImpUpdater<3>;
  template class PoissonBracketImpUpdater<4>;
  template class PoissonBracketImpUpdater<5>;
}

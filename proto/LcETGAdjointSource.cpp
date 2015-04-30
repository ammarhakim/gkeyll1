/**
 * @file	LcETGAdjointSource.cpp
 *
 * @brief	Updater to implement adjoint gyrokinetic equation
 * for the ETG problem
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMathLib.h>
#include <LcETGAdjointSource.h>

namespace Lucee
{
  static const unsigned UPWIND = 0;
  static const unsigned CENTRAL = 1;

// set id for module system
  template <> const char *ETGAdjointSource<4>::id = "ETGAdjointSource4D";
  template <> const char *ETGAdjointSource<5>::id = "ETGAdjointSource5D";

  template <unsigned NDIM>
  ETGAdjointSource<NDIM>::ETGAdjointSource()
    : UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void 
  ETGAdjointSource<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("ETGAdjointSource::readInput: Must specify element to use using 'basis'");

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
        Lucee::Except lce("ETGAdjointSource::readInput: 'fluxType' ");
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
      throw Lucee::Except("ETGAdjointSource::readInput: A jacobianField MUST be supplied!");
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
          throw Lucee::Except("ETGAdjointSource::readInput: updateDirections must be a table with each element < NDIM");
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
          throw Lucee::Except("ETGAdjointSource::readInput: zeroFluxDirections must be a table with each element < NDIM");
      }
    }

    if (tbl.hasNumber("kineticMass"))
      kineticMass = tbl.getNumber("kineticMass");
    else
      throw Lucee::Except(
        "ETGAdjointSource::readInput: Must specify mass using 'kineticMass'");

    if (tbl.hasNumber("eV"))
      eV = tbl.getNumber("eV");
    else
      throw Lucee::Except(
        "ETGAdjointSource::readInput: Must specify units of eV using 'eV'");

    if (tbl.hasNumber("tauOverZi"))
      tauOverZi = tbl.getNumber("tauOverZi");
    else
      throw Lucee::Except(
        "ETGAdjointSource::readInput: Must specify tau/Z_i using 'tauOverZi'");
  }

  template <unsigned NDIM>
  void 
  ETGAdjointSource<NDIM>::initialize()
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

    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex];
    
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
    copyLuceeToEigen(tempVolCoords, volQuad.coordMat);

    // Compute gradients of basis functions evaluated at volume quadrature points
    for (int d = 0; d < updateDirs.size(); d++)
    {
      int dir = updateDirs[d];
      // Each row is a quadrature point; each column is a basis function with derivative applied
      Eigen::MatrixXd derivMatrix = volQuad.interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

      // Store gradients of the same basis function at various quadrature points
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        volQuad.pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);
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
    bigStoredVolMatrices.resize(totalSize);

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

      bigStoredVolMatrices[cellIndex] = massMatrixInv*volQuad.interpMat.transpose();
    }

    massMatrixInv = massMatrix.inverse();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  ETGAdjointSource<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& aCurr = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& stdPotential = this->getInp<Lucee::Field<NDIM, double> >(1);
    const Lucee::Field<NDIM, double>& adjointPotential = this->getInp<Lucee::Field<NDIM, double> >(2);
    const Lucee::Field<NDIM, double>& kineticTemp = this->getInp<Lucee::Field<NDIM, double> >(3);
    const Lucee::Field<NDIM, double>& backgroundF = this->getInp<Lucee::Field<NDIM, double> >(4);
    const Lucee::Field<NDIM, double>& bField4d = this->getInp<Lucee::Field<NDIM, double> >(5);
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

    int nlocal = nodalBasis->getNumNodes();
    int nVolQuad = nodalBasis->getNumGaussNodes();
    int nSurfQuad = nodalBasis->getNumSurfGaussNodes();

    Lucee::ConstFieldPtr<double> aCurrPtr = aCurr.createConstPtr();
    Lucee::ConstFieldPtr<double> stdPotentialPtr = stdPotential.createConstPtr();
    Lucee::ConstFieldPtr<double> adjointPotentialPtr = adjointPotential.createConstPtr();
    Lucee::ConstFieldPtr<double> kineticTempPtr = kineticTemp.createConstPtr();
    Lucee::ConstFieldPtr<double> backgroundFPtr = backgroundF.createConstPtr();
    Lucee::ConstFieldPtr<double> bField4dPtr = bField4d.createConstPtr();
    Lucee::FieldPtr<double> aNewPtr = aNew.createPtr();

    aNew = 0.0; // use aNew to store increment initially

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    double xc[NDIM];

    // Contributions from volume integrals
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      aCurr.setPtr(aCurrPtr, idx);
      aNew.setPtr(aNewPtr, idx);
      stdPotential.setPtr(stdPotentialPtr, idx);
      adjointPotential.setPtr(adjointPotentialPtr, idx);
      kineticTemp.setPtr(kineticTempPtr, idx);
      backgroundF.setPtr(backgroundFPtr, idx);
      bField4d.setPtr(bField4dPtr, idx);

      int cellIndex = volIdxr.getIndex(idx);

      // Compute gradient of adjoint and standard potentials and temperature
      // evaluated at volume quadrature points
      // Each row is a different component/direction
      Eigen::MatrixXd stdPotentialDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nVolQuad);
      Eigen::MatrixXd adjointPotentialDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nVolQuad);
      Eigen::MatrixXd kineticTempDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nVolQuad);
      for (int i = 0; i < nlocal; i++)
      {
        stdPotentialDerivAtQuad += stdPotentialPtr[i]*volQuad.pDiffMatrix[i];
        adjointPotentialDerivAtQuad += adjointPotentialPtr[i]*volQuad.pDiffMatrix[i];
        kineticTempDerivAtQuad += kineticTempPtr[i]*volQuad.pDiffMatrix[i];
      }

      // Interpolate various fields to quadrature points
      Eigen::VectorXd kineticTempVec(nlocal);
      Eigen::VectorXd backgroundFVec(nlocal);
      Eigen::VectorXd bField4dVec(nlocal);
      for (int i = 0; i < nlocal; i++)
      {
        kineticTempVec(i) = kineticTempPtr[i];
        backgroundFVec(i) = backgroundFPtr[i];
        bField4dVec(i) = bField4dPtr[i];
      }
      Eigen::VectorXd kineticTempAtQuad = volQuad.interpMat*kineticTempVec;
      Eigen::VectorXd backgroundFAtQuad = volQuad.interpMat*backgroundFVec;
      Eigen::VectorXd bField4dAtQuad = volQuad.interpMat*bField4dVec;

      grid.setIndex(idx);
      grid.getCentroid(xc);

      //Eigen::VectorXd outputVec(nlocal);

      for (int i = 0; i < nlocal; i++)
      {
        //outputVec(i) = 0.0;
        for (int qp = 0; qp < nVolQuad; qp++)
        {
          // Get velocity coordinate of quadrature point
          double vCoord = xc[2] + volQuad.coordMat(qp,2)*grid.getDx(2)/2.0;
          double muCoord = xc[3] + volQuad.coordMat(qp,3)*grid.getDx(3)/2.0;
          //std::cout << "v1 = " << kineticMass*vCoord*vCoord << std::endl;
          //std::cout << "v2 = " << 2.0*muCoord*bField4dAtQuad(qp) << std::endl << std::endl;

          Eigen::Vector3d bHat(0, 0, 1);
          Eigen::Vector3d gradT(kineticTempDerivAtQuad(0,qp), kineticTempDerivAtQuad(1,qp), 0);
          Eigen::Vector3d gradStdPotential(stdPotentialDerivAtQuad(0,qp),stdPotentialDerivAtQuad(1,qp),0);
          Eigen::Vector3d gradAdjointPotential(adjointPotentialDerivAtQuad(0,qp),adjointPotentialDerivAtQuad(1,qp),0);

          double intermediateResult1 = bHat.cross(gradStdPotential).dot(gradT);
          double intermediateResult2 = bHat.cross(gradAdjointPotential).dot(gradT);

          aNewPtr[i] += volQuad.weights(qp)*bigStoredVolMatrices[cellIndex](i,qp)*
            backgroundFAtQuad(qp)*(intermediateResult1*((kineticMass*vCoord*vCoord + 2.0*muCoord*bField4dAtQuad(qp))/(2.0*kineticTempAtQuad(qp)*eV)
                  -3.0/2.0 + (3.0/2.0)/(1 + tauOverZi))
              - intermediateResult2*kineticMass/(2.0*kineticTempAtQuad(qp)*eV*(1 + tauOverZi)))/kineticTempAtQuad(qp);

          //outputVec(i) += volQuad.weights(qp)*volQuad.interpMat(qp, i)*
          //  backgroundFAtQuad(qp)/bField4dAtQuad(qp)*(intermediateResult1*((kineticMass*vCoord*vCoord + 2.0*muCoord*bField4dAtQuad(qp))/(2.0*kineticTempAtQuad(qp)*eV)
          //        -3.0/2.0 + (3.0/2.0)/(1 + tauOverZi))
          //    - intermediateResult2*kineticMass/(2.0*kineticTempAtQuad(qp)*eV*(1 + tauOverZi)))/kineticTempAtQuad(qp);
        }
      }

      //Eigen::VectorXd realOutputVec = massMatrixInv*outputVec;
      //for (int i = 0; i < nlocal; i++)
      //  aNewPtr[i] = realOutputVec(i);
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

    return Lucee::UpdaterStatus(true, dt);
  }
  
  template <unsigned NDIM>
  void
  ETGAdjointSource<NDIM>::declareTypes()
  {
    // takes two inputs (aCurr, stdPotential, adjointPotential, kineticTemp, backgroundF, bField4d)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns one output, (aNew)
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  ETGAdjointSource<NDIM>::computeNumericalFlux(const Eigen::VectorXd& alphaDotN,
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
    }
  }

  template <unsigned NDIM>
  void
  ETGAdjointSource<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class ETGAdjointSource<4>;
  template class ETGAdjointSource<5>;
}

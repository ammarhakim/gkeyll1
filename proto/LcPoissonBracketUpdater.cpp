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

    if (tbl.hasObject<Lucee::PoissonBracketEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::PoissonBracketEquation>("equation");
    else
    {
      Lucee::Except lce("PoissonBracketUpdater::readInput: Must specify an equation to solve!");
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

    // Check to see if a jacobian field has been supplied
    hasJacobian = false;
    if (tbl.hasBool("hasJacobian"))
      hasJacobian = tbl.getBool("hasJacobian");

    if (hasJacobian == true)
      jacobianField = &tbl.getObject<Lucee::Field<1, double> >("jacobianField");
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
    std::vector<Eigen::MatrixXd> gradStiffnessMatrix(NDIM);
    for (int dir = 0; dir < NDIM; dir++)
    {
      Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, tempMatrix);
      gradStiffnessMatrix[dir] = Eigen::MatrixXd(nlocal, nlocal);

      copyLuceeToEigen(tempMatrix, gradStiffnessMatrix[dir]);
    }

    // get data needed for Gaussian quadrature
    int nVolQuad = nodalBasis->getNumGaussNodes();
    std::vector<double> volWeights(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NDIM);
    volQuad.reset(nVolQuad, nlocal);

    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex];
    
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
    copyLuceeToEigen(tempVolCoords, volQuad.coordMat);

    // Compute gradients of basis functions evaluated at volume quadrature points
    for (int dir = 0; dir < NDIM; dir++)
    {
      // Each row is a quadrature point; each column is a basis function with derivative applied
      Eigen::MatrixXd derivMatrix = volQuad.interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

      // Copy derivatives to different kind type of matrix
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        volQuad.pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);
    }

    // Get data for surface quadrature
    int nSurfQuad = nodalBasis->getNumSurfGaussNodes();
    
    for (int dir = 0; dir < NDIM; dir++)
    {
      // temporary variables
      std::vector<double> tempSurfWeights(nSurfQuad);
      Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal);
      Lucee::Matrix<double> tempSurfCoords(nSurfQuad, NDIM);

      surfLowerQuad[dir].reset(nSurfQuad, nlocal);
      surfUpperQuad[dir].reset(nSurfQuad, nlocal);
      
      // lower surface data
      nodalBasis->getSurfLowerGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfLowerQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex];
      copyLuceeToEigen(tempSurfQuad, surfLowerQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfLowerQuad[dir].coordMat);

      // upper surface data
      nodalBasis->getSurfUpperGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfUpperQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex];
      copyLuceeToEigen(tempSurfQuad, surfUpperQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfUpperQuad[dir].coordMat);
    }

    // Compute gradients of basis functions evaluated at surface quadrature points
    for (int dir = 0; dir < NDIM; dir++)
    {
      // Each row is a quadrature point; each column is a basis function with derivative applied
      Eigen::MatrixXd derivMatrix = surfUpperQuad[dir].interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

      // Use derivMatrix to create matrices for gradients of basis functions
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        surfUpperQuad[dir].pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);

      // Do the same for the lower surface
      derivMatrix = surfLowerQuad[dir].interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
        surfLowerQuad[dir].pDiffMatrix[basisIndex].row(dir) = derivMatrix.col(basisIndex);
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  PoissonBracketUpdater<NDIM>::update(double t)
  {
    std::cout << "update" << std::endl;
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& aCurr = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& hamil = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& aNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    double dt = t-this->getCurrTime();
    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

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
    // Testing
    //if (hasJacobian == true)
    //  Lucee::ConstFieldPtr<double> jacobianPtr = jacobianField->createConstPtr();

    aNew = 0.0; // use aNew to store increment initially

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    // Contributions from volume integrals
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      aCurr.setPtr(aCurrPtr, idx);
      hamil.setPtr(hamilPtr, idx);
      aNew.setPtr(aNewPtr, idx);

      // Compute gradient of hamiltonian
      Eigen::MatrixXd hamilDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nVolQuad);
      for (int i = 0; i < nlocal; i++)
        hamilDerivAtQuad += hamilPtr[i]*volQuad.pDiffMatrix[i];

      Eigen::MatrixXd alpha(NDIM, nVolQuad);
      // TODO: get alpha from appropriate function
      equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, volQuad.interpMat, idx, alpha);

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
      resultVector = massMatrixInv*resultVector;

      for (int i = 0; i < nlocal; i++)
        aNewPtr[i] += resultVector(i);
    }

    // Contributions from surface integrals
    for (int dir = 0; dir < NDIM; dir++)
    {
      // create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq1D(localRgn.deflate(dir));
      // lower and upper bounds of 1D slice. (We need to make sure that flux
      // is computed for one edge outside domain interior)
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

      int idxr[NDIM];
      int idxl[NDIM];

      std::cout << "NDIM = " << NDIM << std::endl;
      std::cout << "sliceLower = " << sliceLower << std::endl << "sliceUpper = " << sliceUpper << std::endl;

      // loop over each 1D slice
      while (seq1D.step())
      {
        seq1D.fillWithIndex(idxr);
        seq1D.fillWithIndex(idxl);
        // loop over each edge
        for (int sliceIndex = sliceLower; sliceIndex < sliceUpper; sliceIndex++)
        { 
          idxr[dir] = sliceIndex;
          idxl[dir] = sliceIndex-1;

          grid.setIndex(idxl);
          double dxL = grid.getDx(dir);
          grid.setIndex(idxr);
          double dxR = grid.getDx(dir);

          double dtdx = 2*dt/(dxL+dxR);

          aCurr.setPtr(aCurrPtr_r, idxr);
          aCurr.setPtr(aCurrPtr_l, idxl);
          // Hamiltonian is continuous, so use left always
          hamil.setPtr(hamilPtr, idxl);
          // Copy data to Eigen vectors
          Eigen::VectorXd rightData(nlocal);
          Eigen::VectorXd leftData(nlocal);
          for (int i = 0; i < nlocal; i++)
          {
            rightData(i) = aCurrPtr_r[i];
            leftData(i) = aCurrPtr_l[i];
          }

          // Compute gradient of hamiltonian at surface nodes
          Eigen::MatrixXd hamilDerivAtQuad = Eigen::MatrixXd::Zero(NDIM, nSurfQuad);
          for (int i = 0; i < nlocal; i++)
            hamilDerivAtQuad += hamilPtr[i]*surfUpperQuad[dir].pDiffMatrix[i];
          // Compute alpha at edge quadrature nodes (making use of alpha dot n being continuous)
          Eigen::MatrixXd alpha(NDIM, nSurfQuad);
          // TODO: get alpha from appropriate function
          equation->computeAlphaAtQuadNodes(hamilDerivAtQuad, surfUpperQuad[dir].interpMat, idx, alpha);
          //std::cout << "rows = " << surfUpperQuad[dir].interpMat.rows() << std::endl;

          // Construct normal vector
          Eigen::VectorXd normalVec = Eigen::VectorXd::Zero(NDIM);
          normalVec(dir) = 1.0;
          
          // Calculate alphaDotN at all quadrature points
          Eigen::VectorXd alphaDotN = normalVec.transpose()*alpha;

          // Compute numerical flux
          Eigen::VectorXd numericalFluxAtQuad(nSurfQuad);
          computeNumericalFlux(alphaDotN, surfUpperQuad[dir].interpMat*leftData,
              surfLowerQuad[dir].interpMat*rightData, numericalFluxAtQuad);

          // CFL adjustment
          double maxSpeed = alphaDotN.cwiseAbs().maxCoeff();
          if (dtdx*maxSpeed > cfla)
          {
            cfla = dtdx*maxSpeed;
            // Check if time-step was too large and return a suggestion with correct time-step
            if (cfla > cflm)
            {
              std::cout << "break" << std::endl;
              return Lucee::UpdaterStatus(false, dt*cfl/cfla);
            }
          }

          // Compute the surface integrals required and multiply by inverse of mass matrix
          Eigen::VectorXd upperResultVector = massMatrixInv*surfUpperQuad[dir].interpMat.transpose()*
            numericalFluxAtQuad.cwiseProduct(surfUpperQuad[dir].weights);
          Eigen::VectorXd lowerResultVector = massMatrixInv*surfLowerQuad[dir].interpMat.transpose()*
            numericalFluxAtQuad.cwiseProduct(surfLowerQuad[dir].weights);

          // Accumulate to solution
          aNew.setPtr(aNewPtr_l, idxr);
          aNew.setPtr(aNewPtr_r, idxl);
          for (int i = 0; i < nlocal; i++)
          {
            aNewPtr_l[i] -= upperResultVector(i);
            aNewPtr_r[i] += lowerResultVector(i);
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
  PoissonBracketUpdater<NDIM>::declareTypes()
  {
    // takes two inputs (aCurr, hamil)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns one output, (aNew)
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  PoissonBracketUpdater<NDIM>::computeNumericalFlux(const Eigen::VectorXd& alphaDotN,
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

/**
 * @file	LcLenardBernsteinDrag5DUpdater.cpp
 *
 * @brief	Updater to compute the drag term in the L-B collision operator.
 * df/dt = alpha*d[(v-uIn)f]/dv where uIn is a function of (x,y,z) (drift velocity)
 * Used for 3D2V SOL problem.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLenardBernsteinDrag5DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *LenardBernsteinDrag5DUpdater::id = "LenardBernsteinDrag5DUpdater";

  LenardBernsteinDrag5DUpdater::LenardBernsteinDrag5DUpdater()
    : UpdaterIfc()
  {
  }

  void 
  LenardBernsteinDrag5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("LenardBernsteinDrag5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("LenardBernsteinDrag5DUpdater::readInput: Must specify element to use using 'basis3d'");
    
    // cfl number
    cfl = tbl.getNumber("cfl");
    // use slightly large max CFL to avoid thrashing around
    cflm = 1.1*cfl;
    
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    if (tbl.hasFunction("alpha"))
      fnRef = tbl.getFunctionRef("alpha");
    else
      throw Lucee::Except("LenardBernsteinDrag5DUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void 
  LenardBernsteinDrag5DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();
    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<5> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[5];
    seq.fillWithIndex(idx);
    nodalBasis5d->setIndex(idx);
    
    int nlocal5d = nodalBasis5d->getNumNodes();
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    int nSurfQuad5d = nodalBasis5d->getNumSurfGaussNodes();

    // Get mass and grad-stiffness matrices and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal5d, nlocal5d);
    Lucee::Matrix<double> gradStiffMatrixLucee(nlocal5d, nlocal5d);

    Eigen::MatrixXd massMatrix(nlocal5d, nlocal5d);
    nodalBasis5d->getMassMatrix(massMatrixLucee);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    // Element 0 is dir 3, Element 1 is dir 4
    std::vector<Eigen::MatrixXd> gradStiffMatrix(2, Eigen::MatrixXd(nlocal5d, nlocal5d));
    nodalBasis5d->getGradStiffnessMatrix(3, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix[0]);
    nodalBasis5d->getGradStiffnessMatrix(4, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix[1]);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpVolMatrixLucee(nVolQuad5d, nlocal5d);
    Lucee::Matrix<double> interpSurfMatrixLower5dLucee(nSurfQuad5d, nlocal5d);
    Lucee::Matrix<double> interpSurfMatrixUpper5dLucee(nSurfQuad5d, nlocal5d);

    Lucee::Matrix<double> gaussVolOrdinatesLucee(nVolQuad5d, 5);
    Lucee::Matrix<double> gaussSurfOrdinatesLucee(nSurfQuad5d, 5);
    gaussVolWeights5d = std::vector<double>(nVolQuad5d);
    gaussVolOrdinates5d = Eigen::MatrixXd(nVolQuad5d, 5);
    
    // Get the interpolation matrix for the volume quadrature points
    interpVolMatrix5d = Eigen::MatrixXd(nVolQuad5d, nlocal5d);
    nodalBasis5d->getGaussQuadData(interpVolMatrixLucee, gaussVolOrdinatesLucee, gaussVolWeights5d);
    copyLuceeToEigen(interpVolMatrixLucee, interpVolMatrix5d);
    copyLuceeToEigen(gaussVolOrdinatesLucee, gaussVolOrdinates5d);

    interpSurfMatrixLower5d.resize(2, Eigen::MatrixXd(nSurfQuad5d, nlocal5d));
    interpSurfMatrixUpper5d.resize(2, Eigen::MatrixXd(nSurfQuad5d, nlocal5d));
    gaussSurfOrdinates5d.resize(2, Eigen::MatrixXd(nSurfQuad5d, 5));
    gaussSurfWeights5d.resize(2, std::vector<double>(nSurfQuad5d));
    // Get the interpolation matrix for the upper surface quadrature points.
    // Quadrature location and weights will be overwritten but it doesn't matter for these purposes
    nodalBasis5d->getSurfUpperGaussQuadData(3, interpSurfMatrixUpper5dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights5d[0]);
    copyLuceeToEigen(interpSurfMatrixUpper5dLucee, interpSurfMatrixUpper5d[0]);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis5d->getSurfLowerGaussQuadData(3, interpSurfMatrixLower5dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights5d[0]);
    copyLuceeToEigen(interpSurfMatrixLower5dLucee, interpSurfMatrixLower5d[0]);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates5d[0]);

    // Do the same thing as above in the mu direction
    nodalBasis5d->getSurfUpperGaussQuadData(4, interpSurfMatrixUpper5dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights5d[1]);
    copyLuceeToEigen(interpSurfMatrixUpper5dLucee, interpSurfMatrixUpper5d[1]);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis5d->getSurfLowerGaussQuadData(4, interpSurfMatrixLower5dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights5d[1]);
    copyLuceeToEigen(interpSurfMatrixLower5dLucee, interpSurfMatrixLower5d[1]);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates5d[1]);

    // Compute and store inverse of mass matrix
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();
    // Compute derivative of basis function evaluated at volume quadrature points
    basisDerivAtVolQuad.resize(2);
    basisDerivAtVolQuad[0] = interpVolMatrix5d*massMatrixInv*gradStiffMatrix[0].transpose();
    basisDerivAtVolQuad[1] = interpVolMatrix5d*massMatrixInv*gradStiffMatrix[1].transpose();
    // Take transpose so same basis function evaluated at diff quad points in the same row
    // Need to use transposeInPlace because otherwise there will be a bug!
    basisDerivAtVolQuad[0].transposeInPlace();
    basisDerivAtVolQuad[1].transposeInPlace();

    // Pre-multiply by inverse mass matrix
    basisDerivAtVolQuad[0] = massMatrixInv*basisDerivAtVolQuad[0];
    basisDerivAtVolQuad[1] = massMatrixInv*basisDerivAtVolQuad[1];

    // Precompute two matrices for surface integrals
    surfIntegralMatrixLower5d.resize(2);
    surfIntegralMatrixUpper5d.resize(2);
    // Pre-multiply by inverse mass matrix
    surfIntegralMatrixLower5d[0] = massMatrixInv*interpSurfMatrixLower5d[0].transpose();
    surfIntegralMatrixUpper5d[0] = massMatrixInv*interpSurfMatrixUpper5d[0].transpose();
    surfIntegralMatrixLower5d[1] = massMatrixInv*interpSurfMatrixLower5d[1].transpose();
    surfIntegralMatrixUpper5d[1] = massMatrixInv*interpSurfMatrixUpper5d[1].transpose();

    nodalBasis3d->setIndex(idx[0], idx[1], idx[2]);
    int nlocal3d = nodalBasis3d->getNumNodes();
    // Compute mom0Vector for cell-average calculation
    Lucee::Matrix<double> massMatrix3dLucee(nlocal3d, nlocal3d);
    nodalBasis3d->getMassMatrix(massMatrix3dLucee);
    Eigen::MatrixXd massMatrix3d(nlocal3d, nlocal3d);
    copyLuceeToEigen(massMatrix3dLucee, massMatrix3d);
    mom0Vector = massMatrix3d.colwise().sum();
    // Get the 3d interpolation matrix
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> gaussVolWeights3d(nVolQuad3d);
    Lucee::Matrix<double> gaussVolOrdinates3dLucee(nVolQuad3d, 3);
    Lucee::Matrix<double> interpVolMatrix3dLucee(nVolQuad3d, nlocal3d);
    interpVolMatrix3d = Eigen::MatrixXd(nVolQuad3d, nlocal3d);
    nodalBasis3d->getGaussQuadData(interpVolMatrix3dLucee, gaussVolOrdinates3dLucee, gaussVolWeights3d);
    copyLuceeToEigen(interpVolMatrix3dLucee, interpVolMatrix3d);
  }

  Lucee::UpdaterStatus 
  LenardBernsteinDrag5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& fIn = this->getInp<Lucee::Field<5, double> >(0);
    // Drift velocity uIn(x).
    const Lucee::Field<3, double>& uIn = this->getInp<Lucee::Field<3, double> >(1);
    // Number density at node
    const Lucee::Field<3, double>& nIn = this->getInp<Lucee::Field<3, double> >(2);
    // Temperature in joules
    const Lucee::Field<3, double>& temperatureIn = this->getInp<Lucee::Field<3, double> >(3);
    // Dimensionally correct number density from weighted moment calculation
    const Lucee::Field<3, double>& numDensityIn = this->getInp<Lucee::Field<3, double> >(4);
    // Output distribution function
    Lucee::Field<5, double>& fOut = this->getOut<Lucee::Field<5, double> >(0);

    int nlocal5d = nodalBasis5d->getNumNodes();
    int nlocal3d = nodalBasis3d->getNumNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();

    double dt = t-this->getCurrTime();

    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();

    double cfla = 0.0; // maximum CFL number

    Lucee::ConstFieldPtr<double> fInPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInPtrLeft = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> uInPtr = uIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nInPtr = nIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureInPtr = temperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> numDensityInPtr = numDensityIn.createConstPtr();
    Lucee::FieldPtr<double> fOutPtr = fOut.createPtr();
    Lucee::FieldPtr<double> fOutPtrLeft = fOut.createPtr();

    fOut = 0.0; // use fOut to store increment initially
    int idx[5];
    double cellCentroid[5];
    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    // Get alpha. Need to scale by n/T^(3/2) in this class
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    alpha = resultVector[0];

    // Volume integral contribution for both v-parallel and mu directions
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);
      
      fIn.setPtr(fInPtr, idx);
      uIn.setPtr(uInPtr, idx[0], idx[1], idx[2]);
      nIn.setPtr(nInPtr, idx[0], idx[1], idx[2]);
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);
      fOut.setPtr(fOutPtr, idx);
      
      Eigen::VectorXd uVals(nlocal3d);
      Eigen::VectorXd nVals(nlocal3d);
      // Copy values of uIn into uVals vector
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        uVals(nodeIndex) = uInPtr[nodeIndex];
        nVals(nodeIndex) = nInPtr[nodeIndex];
      }
      // Interpolate uIn to 3D quadrature points 
      Eigen::VectorXd uAtQuad = interpVolMatrix3d*uVals;
      Eigen::VectorXd nAtQuad = interpVolMatrix3d*nVals;
       
      // Compute distribution function at quadrature points
      Eigen::VectorXd fVals(nlocal5d);
      for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        fVals(nodeIndex) = fInPtr[nodeIndex];
      
      Eigen::VectorXd fVolQuad = interpVolMatrix5d*fVals;

      Eigen::VectorXd paraQuad(nVolQuad5d);
      Eigen::VectorXd muQuad(nVolQuad5d);
      // Compute integrands for gaussian quadrature excluding basis function derivatives
      for (int quadIndex = 0; quadIndex < nVolQuad5d; quadIndex++)
      {
        double paraCoord = cellCentroid[3] + 0.5*grid.getDx(3)*gaussVolOrdinates5d(quadIndex,3);
        double muCoord = cellCentroid[4] + 0.5*grid.getDx(4)*gaussVolOrdinates5d(quadIndex,4);
        paraQuad(quadIndex) = gaussVolWeights5d[quadIndex]*fVolQuad(quadIndex)*
          (paraCoord - uAtQuad(quadIndex % nVolQuad3d)/nAtQuad(quadIndex % nVolQuad3d));
        muQuad(quadIndex) = gaussVolWeights5d[quadIndex]*fVolQuad(quadIndex)*2*muCoord;
      }

      // Evaluate integral using gaussian quadrature (represented as matrix-vector multiply)
      Eigen::VectorXd volIntegralResult = (basisDerivAtVolQuad[0]*paraQuad +
        basisDerivAtVolQuad[1]*muQuad);
      
      for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        fOutPtr[nodeIndex] -= volIntegralResult(nodeIndex);
    }

    // Contributions from surface integrals
    // Create sequencer to loop over *each* 4D slice in '3' direction
    Lucee::RowMajorSequencer<5> seqLowerDim(localRgn.deflate(3));

    // To set indices on left and right cells of a boundary
    int idxr[5];
    int idxl[5];
    
    int ivLower = localRgn.getLower(3);
    // Need one edge outside domain interior
    int ivUpper = localRgn.getUpper(3)+1;

    // Zero flux BC's. Make sure we skip surface integrals on
    // surfaces that lie on a zero-flux boundary
    if (ivLower == globalRgn.getLower(3))
      ivLower = globalRgn.getLower(3)+1;
    if (ivUpper == globalRgn.getUpper(3)+1)
      ivUpper = globalRgn.getUpper(3);

    // Flag to keep track of whether or not CFL limit was exceeded.
    // Don't exit immediately upon false because a more strict limit might appear
    bool needToRetakeStep = false;
    
    // Loop over each 4D cell in (x,y,z,mu) space
    while (seqLowerDim.step())
    {
      seqLowerDim.fillWithIndex(idxr);
      seqLowerDim.fillWithIndex(idxl);

      // Compute alpha scale factor n/Te^(3/2)
      temperatureIn.setPtr(temperatureInPtr, idxr[0], idxr[1], idxr[2]);
      numDensityIn.setPtr(numDensityInPtr, idxr[0], idxr[1], idxr[2]);
      Eigen::VectorXd temperatureVec(nlocal3d);
      Eigen::VectorXd numDensityVec(nlocal3d);
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        temperatureVec(nodeIndex) = temperatureInPtr[nodeIndex];
        numDensityVec(nodeIndex) = numDensityInPtr[nodeIndex];
      }
      grid.setIndex(idxr);
      double volume = grid.getDx(0)*grid.getDx(1)*grid.getDx(2);
      double averageTemperature = mom0Vector.dot(temperatureVec)/volume;
      double averageN = mom0Vector.dot(numDensityVec)/volume;

      // Compute uIn at volume quadrature points (again)
      uIn.setPtr(uInPtr, idxr[0], idxr[1], idxr[2]);
      nIn.setPtr(nInPtr, idxr[0], idxr[1], idxr[2]);
      Eigen::VectorXd uVals(nlocal3d);
      Eigen::VectorXd nVals(nlocal3d);
      // Copy values of uIn into uVals vector
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        uVals(nodeIndex) = uInPtr[nodeIndex];
        nVals(nodeIndex) = nInPtr[nodeIndex];
      }
      // Interpolate uIn to 3D quadrature points 
      Eigen::VectorXd uAtQuad = interpVolMatrix3d*uVals;
      Eigen::VectorXd nAtQuad = interpVolMatrix3d*nVals;

      // Loop over each edge in vParallel
      for (int i = ivLower; i < ivUpper; i++)
      {
        idxr[3] = i; // cell right of edge
        idxl[3] = i-1; // cell left of edge
        
        fIn.setPtr(fInPtr, idxr);
        fIn.setPtr(fInPtrLeft, idxl);
        
        grid.setIndex(idxl);
        double dxL = grid.getDx(3);
        grid.setIndex(idxr);
        double dxR = grid.getDx(3);

        // Copy fIn data to Eigen vectors for matrix multiplications
        Eigen::VectorXd fLeft(nlocal5d);
        Eigen::VectorXd fRight(nlocal5d);
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        {
          fLeft(nodeIndex) = fInPtrLeft[nodeIndex];
          fRight(nodeIndex) = fInPtr[nodeIndex];
        }

        // Evaluate fLeft and fRight at quadrature nodes on surface
        Eigen::VectorXd fLeftSurfEvals = interpSurfMatrixUpper5d[0]*fLeft;
        Eigen::VectorXd fRightSurfEvals = interpSurfMatrixLower5d[0]*fRight;

        // Loop over quadrature points on the edge
        Eigen::VectorXd surfIntegralFluxes(gaussSurfWeights5d[0].size());

        // Need to know center of right cell to figure out the global velocity coordinate
        grid.getCentroid(cellCentroid);
        double paraCoord = cellCentroid[3] - 0.5*grid.getDx(3);
        
        for (int quadIndex = 0; quadIndex < gaussSurfWeights5d[0].size(); quadIndex++)
        {
          // Compute Lax flux at each surface quadrature point
          double numFlux;
          if (paraCoord - uAtQuad(quadIndex % nVolQuad3d)/nAtQuad(quadIndex % nVolQuad3d) > 0)
            numFlux = (paraCoord - uAtQuad(quadIndex % nVolQuad3d)/nAtQuad(quadIndex % nVolQuad3d))*
              fRightSurfEvals(quadIndex);
          else numFlux = (paraCoord - uAtQuad(quadIndex % nVolQuad3d)/nAtQuad(quadIndex % nVolQuad3d))*
            fLeftSurfEvals(quadIndex);

          // Keep track of max CFL number
          cfla = std::max(cfla, std::abs(alpha*averageN/(averageTemperature*sqrt(averageTemperature))*
                (paraCoord-uAtQuad(quadIndex % nVolQuad3d)/nAtQuad(quadIndex % nVolQuad3d))*dt/grid.getDx(3)));
          //cfla = std::max(cfla, std::abs(alpha*(paraCoord-uAtQuad(quadIndex % nVolQuad3d))*dt/(0.5*(dxL+dxR))));
          // Time-step was too large: return a suggestion with correct time-step
          if (cfla > cflm)
            needToRetakeStep = true;

          // Store result of weight*(v-uIn)*f at this location in a vector
          surfIntegralFluxes(quadIndex) = gaussSurfWeights5d[0][quadIndex]*numFlux;
        }

        // Compute all surface integrals using a matrix multiply.
        // Each row in result is a basis function times flux projection
        // Then inverse of mass matrix is multiplied to find appropriate increments
        Eigen::VectorXd leftSurfIntegralResult = surfIntegralMatrixUpper5d[0]*surfIntegralFluxes;
        Eigen::VectorXd rightSurfIntegralResult = surfIntegralMatrixLower5d[0]*surfIntegralFluxes;

        // Update left cell connected to edge with flux on face
        fOut.setPtr(fOutPtrLeft, idxl);
        // Update right cell connected to edge with flux on face
        fOut.setPtr(fOutPtr, idxr);

        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        {
          // Plus sign for left cell since outward normal is in +v direction
          fOutPtrLeft[nodeIndex] += leftSurfIntegralResult(nodeIndex);
          // Minus sign for right cell since outward normal is in -v direction
          fOutPtr[nodeIndex]  -= rightSurfIntegralResult(nodeIndex);
        }
      }
    }

    // Contributions from surface integrals
    // Create sequencer to loop over *each* 4D slice in '4' direction
    seqLowerDim = Lucee::RowMajorSequencer<5>(localRgn.deflate(4));

    int iMuLower = localRgn.getLower(4);
    // Need one edge outside domain interior
    int iMuUpper = localRgn.getUpper(4)+1;

    // Zero flux BC's
    if (iMuLower == globalRgn.getLower(4))
      iMuLower = globalRgn.getLower(4)+1;
    if (iMuUpper == globalRgn.getUpper(4)+1)
      iMuUpper = globalRgn.getUpper(4);
  
    // Loop over each 4D slice in (x,y,z,v)
    while (seqLowerDim.step())
    {
      seqLowerDim.fillWithIndex(idxr);
      seqLowerDim.fillWithIndex(idxl);

      // Compute alpha scale factor n/Te^(3/2)
      temperatureIn.setPtr(temperatureInPtr, idxr[0], idxr[1], idxr[2]);
      numDensityIn.setPtr(numDensityInPtr, idxr[0], idxr[1], idxr[2]);
      Eigen::VectorXd temperatureVec(nlocal3d);
      Eigen::VectorXd numDensityVec(nlocal3d);
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        temperatureVec(nodeIndex) = temperatureInPtr[nodeIndex];
        numDensityVec(nodeIndex) = numDensityInPtr[nodeIndex];
      }
      grid.setIndex(idxr);
      double volume = grid.getDx(0)*grid.getDx(1)*grid.getDx(2);
      double averageTemperature = mom0Vector.dot(temperatureVec)/volume;
      double averageN = mom0Vector.dot(numDensityVec)/volume;

      // Loop over each edge in mu
      for (int i = iMuLower; i < iMuUpper; i++)
      {
        idxr[4] = i; // cell right of edge
        idxl[4] = i-1; // cell left of edge
        
        fIn.setPtr(fInPtr, idxr);
        fIn.setPtr(fInPtrLeft, idxl);
        
        grid.setIndex(idxl);
        double dxL = grid.getDx(4);
        grid.setIndex(idxr);
        double dxR = grid.getDx(4);

        // Copy fIn data to Eigen vectors for matrix multiplications
        Eigen::VectorXd fLeft(nlocal5d);
        Eigen::VectorXd fRight(nlocal5d);
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        {
          fLeft(nodeIndex) = fInPtrLeft[nodeIndex];
          fRight(nodeIndex) = fInPtr[nodeIndex];
        }

        // Evaluate fLeft and fRight at quadrature nodes on surface
        Eigen::VectorXd fLeftSurfEvals = interpSurfMatrixUpper5d[1]*fLeft;
        Eigen::VectorXd fRightSurfEvals = interpSurfMatrixLower5d[1]*fRight;

        // Loop over quadrature points on the edge
        Eigen::VectorXd surfIntegralFluxes(gaussSurfWeights5d[1].size());

        // Need to know center of right cell to figure out the global velocity coordinate
        grid.getCentroid(cellCentroid);
        double muCoord = cellCentroid[4] - 0.5*grid.getDx(4);
        
        // Keep track of max CFL number
        cfla = std::max(cfla, std::abs(alpha*averageN/(averageTemperature*sqrt(averageTemperature))*
          2*muCoord*dt/grid.getDx(4)));
        // Time-step was too large: return a suggestion with correct time-step
        if (cfla > cflm)
          needToRetakeStep = true;

        for (int quadIndex = 0; quadIndex < gaussSurfWeights5d[1].size(); quadIndex++)
        {
          // Compute upwind flux at each surface quadrature point (mu will always be positive)
          double numFlux = 2*muCoord*fRightSurfEvals(quadIndex);

          // Store result of weight*mu*f at this location in a vector
          surfIntegralFluxes(quadIndex) = gaussSurfWeights5d[1][quadIndex]*numFlux;
        }

        // Compute all surface integrals using a matrix multiply.
        // Each row in result is a basis function times flux projection
        // Then inverse of mass matrix is multiplied to find appropriate increments
        Eigen::VectorXd leftSurfIntegralResult = surfIntegralMatrixUpper5d[1]*surfIntegralFluxes;
        Eigen::VectorXd rightSurfIntegralResult = surfIntegralMatrixLower5d[1]*surfIntegralFluxes;

        // Update left cell connected to edge with flux on face
        fOut.setPtr(fOutPtrLeft, idxl);
        // Update right cell connected to edge with flux on face
        fOut.setPtr(fOutPtr, idxr);

        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        {
          // Plus sign for left cell since outward normal is in +v direction
          fOutPtrLeft[nodeIndex] += leftSurfIntegralResult(nodeIndex);
          // Minus sign for right cell since outward normal is in -v direction
          fOutPtr[nodeIndex] -= rightSurfIntegralResult(nodeIndex);
        }
      }
    }

    if (needToRetakeStep == true)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

    
    seq.reset();
    // Final sweep, update solution with forward Euler step
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fOut.setPtr(fOutPtr, idx);
      // Compute alpha scale factor n/Te^(3/2)
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);
      Eigen::VectorXd temperatureVec(nlocal3d);
      Eigen::VectorXd numDensityVec(nlocal3d);
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        temperatureVec(nodeIndex) = temperatureInPtr[nodeIndex];
        numDensityVec(nodeIndex) = numDensityInPtr[nodeIndex];
      }
      grid.setIndex(idx);
      double volume = grid.getDx(0)*grid.getDx(1)*grid.getDx(2);
      double averageTemperature = mom0Vector.dot(temperatureVec)/volume;
      double averageN = mom0Vector.dot(numDensityVec)/volume;
      
      if (onlyIncrement == false)
      {
        fIn.setPtr(fInPtr, idx);
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          fOutPtr[nodeIndex] = fInPtr[nodeIndex] + dt*alpha*averageN/(averageTemperature*sqrt(averageTemperature))*
            fOutPtr[nodeIndex];
      }
      else
      {
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          fOutPtr[nodeIndex] = alpha*averageN/(averageTemperature*sqrt(averageTemperature))*
            fOutPtr[nodeIndex];
      }
    }
    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDrag5DUpdater::declareTypes()
  {
    // Input: Distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: <v> mean velocity, need to divide by <1>
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: <1> number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: dimensionally correct number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // returns one output (fNew)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  LenardBernsteinDrag5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
    {
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
      {
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
      }
    }
  }

  void
  LenardBernsteinDrag5DUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("LenardBernsteinDrag5DUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'alpha' "
          << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      lce << "[" << err << "]";
      throw lce;
    }
    // fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("LenardBernsteinDrag5DUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}

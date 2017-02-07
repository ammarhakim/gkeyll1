/**
 * @file	LcLenardBernsteinDiffAlternate3D2VUpdater.cpp
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator for 3D2V problems.
 * This updater is different from LcLenardBernsteinDiff5DUpdater in its use of a 2D recovery method
 * to evaluate diffusion terms
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLenardBernsteinDiffAlternate3D2VUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

//etc includes
#include <quadrule.hpp>

namespace Lucee
{
  using namespace Eigen;
  const char *LenardBernsteinDiffAlternate3D2VUpdater::id = "LenardBernsteinDiffAlternate3D2VUpdater";

  LenardBernsteinDiffAlternate3D2VUpdater::LenardBernsteinDiffAlternate3D2VUpdater()
  {
  }

  LenardBernsteinDiffAlternate3D2VUpdater::~LenardBernsteinDiffAlternate3D2VUpdater()
  {
  }

  void
  LenardBernsteinDiffAlternate3D2VUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::readInput: Must specify element to use using 'basis5d'");
 
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::readInput: Must specify element to use using 'basis3d'");
 

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::readInput: Must specify element to use using 'basis2d'");
 
    // CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::readInput: Must specify speciesMass");

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::readInput: Must specify polyOrder");

    if (tbl.hasFunction("alpha"))
      fnRef = tbl.getFunctionRef("alpha");
    else
      throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void
  LenardBernsteinDiffAlternate3D2VUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();
    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<5> seq5d(localRgn);
    seq5d.step(); // just to get to first index
    int idx[5];
    seq5d.fillWithIndex(idx);
    nodalBasis5d->setIndex(idx);
    grid.setIndex(idx);
    nodalBasis2d->setIndex(idx[0],idx[1]);

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    // Get a copy of the nodal coordinates
    Lucee::Matrix<double> nodeCoordsLucee(nlocal5d, 5);
    nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
    Eigen::MatrixXd nodeCoords(nlocal5d, 5);
    copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

    double dxMin = grid.getDx(0);
    for (int d = 1; d < 3; d++)
      dxMin = std::min(dxMin, grid.getDx(d));

    // Find all nodes that share the same location as node zero
    // Will eventually need to do this at all nodes on lower surface
    // but can get away with doing this at one point for linear element
    nodalStencil = std::vector<int>(nlocal5d);
    int stencilIndex = 0;
    for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
    {
      if (sameConfigCoords(0, nodeIndex, dxMin, nodeCoords) == true)
      {
        nodalStencil[stencilIndex] = nodeIndex;
        stencilIndex++;
      }
    }
    nodalStencil.resize(stencilIndex);

    // Compute inverse matrix for premultiplication
    Lucee::Matrix<double> massMatrixLucee(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(massMatrixLucee);
    Eigen::MatrixXd massMatrix(nlocal2d, nlocal2d);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    // Compute mass matrix inverse (remember to remove computational grid scaling)
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();

    int nSurfQuad = nodalBasis2d->getNumSurfGaussNodes();
    int nVolQuad = nodalBasis2d->getNumGaussNodes();
    // Allocate temporary matrices
    std::vector<double> volWeights(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal2d);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NC2);
    Lucee::Matrix<double> tempGradStiffness(nlocal2d, nlocal2d);
    Eigen::MatrixXd gradStiffness(nlocal2d, nlocal2d);
    // Initialize volume quadrature struct
    volQuad.reset(nVolQuad, nlocal2d, NC2);
    // Get volume quadrature matrix
    nodalBasis2d->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
    copyLuceeToEigen(tempVolCoords, volQuad.coordMat);

    surfLowerUpdateMatrix.resize(2);
    surfLowerUpdateMatrixDeriv.resize(2);
    surfUpperUpdateMatrix.resize(2);
    surfUpperUpdateMatrixDeriv.resize(2);
    volumeUpdateMatrix.resize(2);
    // Get data for surface quadrature
    for (int dir = 0; dir < 2; dir++)
    {
      // temporary variables
      std::vector<double> tempSurfWeights(nSurfQuad);
      Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal2d);
      Lucee::Matrix<double> tempSurfCoords(nSurfQuad, NC2);

      // Reset surface quadrature structures
      surfLowerQuad[dir].reset(nSurfQuad, nlocal2d, NC2);
      surfUpperQuad[dir].reset(nSurfQuad, nlocal2d, NC2);
      
      // Compute scale factor to undo scaling of weights
      double surfAreaComputational = 1.0;
      for (int dirIndex = 0; dirIndex < 2; dirIndex++)
      {
        if (dirIndex != dir)
          surfAreaComputational = surfAreaComputational*grid.getDx(dirIndex);
      }

      // Get grad stiffness matrix for this direction
      nodalBasis2d->getGradStiffnessMatrix(dir, tempGradStiffness);
      copyLuceeToEigen(tempGradStiffness, gradStiffness);

      // Compute volume update matrix (really just needed for mu. will need to be modified for v)
      // Need to scale by 1/dv in update
      volumeUpdateMatrix[dir] = grid.getDx(0)*grid.getDx(1)/surfAreaComputational*massMatrixInv*gradStiffness;

      // lower surface data
      nodalBasis2d->getSurfLowerGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfLowerQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex]/surfAreaComputational;
      copyLuceeToEigen(tempSurfQuad, surfLowerQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfLowerQuad[dir].coordMat);

      // each column is a basis function with derivative applied
      // grid.getDx(dir) multiplied on to remove grid scale. Proper value needs to be multiplied in update()
      Eigen::MatrixXd derivMatrixLower = grid.getDx(dir)*surfLowerQuad[dir].interpMat*massMatrixInv*gradStiffness.transpose();

      // Compute lower update matrix
      surfLowerUpdateMatrix[dir] = grid.getDx(0)*grid.getDx(1)*massMatrixInv*surfLowerQuad[dir].interpMat.transpose();
      // Multiply each row of surfLowerUpdateMatrix by the integration weights
      for (int rowIndex = 0; rowIndex < surfLowerUpdateMatrix[dir].rows(); rowIndex++)
        for (int colIndex = 0; colIndex < surfLowerUpdateMatrix[dir].cols(); colIndex++)
          surfLowerUpdateMatrix[dir](rowIndex,colIndex) *= surfLowerQuad[dir].weights(colIndex);

      // Compute lower update matrix derivative
      surfLowerUpdateMatrixDeriv[dir] = grid.getDx(0)*grid.getDx(1)*massMatrixInv*derivMatrixLower.transpose();
      // Multiply each row of surfLowerUpdateMatrix by the integration weights
      for (int rowIndex = 0; rowIndex < surfLowerUpdateMatrixDeriv[dir].rows(); rowIndex++)
        for (int colIndex = 0; colIndex < surfLowerUpdateMatrixDeriv[dir].cols(); colIndex++)
          surfLowerUpdateMatrixDeriv[dir](rowIndex,colIndex) *= surfLowerQuad[dir].weights(colIndex);

      // upper surface data
      nodalBasis2d->getSurfUpperGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfUpperQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex]/surfAreaComputational;
      copyLuceeToEigen(tempSurfQuad, surfUpperQuad[dir].interpMat);
      copyLuceeToEigen(tempSurfCoords, surfUpperQuad[dir].coordMat);
      
      Eigen::MatrixXd derivMatrixUpper = grid.getDx(dir)*surfUpperQuad[dir].interpMat*massMatrixInv*gradStiffness.transpose();

      // Compute upper update matrix
      surfUpperUpdateMatrix[dir] = grid.getDx(0)*grid.getDx(1)*massMatrixInv*surfUpperQuad[dir].interpMat.transpose();
      // Multiply each row of surfUpperUpdateMatrix by the integration weights
      for (int rowIndex = 0; rowIndex < surfUpperUpdateMatrix[dir].rows(); rowIndex++)
        for (int colIndex = 0; colIndex < surfUpperUpdateMatrix[dir].cols(); colIndex++)
          surfUpperUpdateMatrix[dir](rowIndex,colIndex) *= surfUpperQuad[dir].weights(colIndex);

      // Compute upper update matrix derivative
      surfUpperUpdateMatrixDeriv[dir] = grid.getDx(0)*grid.getDx(1)*massMatrixInv*derivMatrixUpper.transpose();
      // Multiply each row of surfUpperUpdateMatrix by the integration weights
      for (int rowIndex = 0; rowIndex < surfUpperUpdateMatrixDeriv[dir].rows(); rowIndex++)
        for (int colIndex = 0; colIndex < surfUpperUpdateMatrixDeriv[dir].cols(); colIndex++)
          surfUpperUpdateMatrixDeriv[dir](rowIndex,colIndex) *= surfUpperQuad[dir].weights(colIndex);
    }

    // Construct recovery quadrature matrices
    // First figure out quadrature points in 1d
    int numGaussPoints1d = polyOrder + 1;
    gaussPoints1d.resize(numGaussPoints1d);
    std::vector<double> gaussWeights1d(numGaussPoints1d);
    legendre_set(numGaussPoints1d, &gaussPoints1d[0], &gaussWeights1d[0]);

    double refCoord[2];
    std::vector<double> basisAtPoint(nlocal2d);
    recoveryInterpMatVPara.resize(nSurfQuad);
    recoveryInterpMatMu.resize(nSurfQuad);
    // At each surface quadrature point, calculate recovery quadrature points
    for (int surfQuadIndex = 0; surfQuadIndex < nSurfQuad; surfQuadIndex++)
    {
      recoveryInterpMatVPara[surfQuadIndex] = Eigen::MatrixXd(numGaussPoints1d,nlocal2d);
      recoveryInterpMatMu[surfQuadIndex] = Eigen::MatrixXd(numGaussPoints1d,nlocal2d);

      for (int quadIndex = 0; quadIndex < numGaussPoints1d; quadIndex++)
      {
        // Fill out coordinate of quadrature point for vPara
        refCoord[0] = gaussPoints1d[quadIndex];
        refCoord[1] = surfUpperQuad[0].coordMat(surfQuadIndex,1);
        // Evaluate all basis functions at point
        nodalBasis2d->evalBasis(refCoord,basisAtPoint);
        // Copy to recoveryInterpMatVPara
        for (int basisIndex = 0; basisIndex < nlocal2d; basisIndex++)
          recoveryInterpMatVPara[surfQuadIndex](quadIndex,basisIndex) = basisAtPoint[basisIndex];

        // Fill out coordinate of quadrature point for mu
        refCoord[0] = surfUpperQuad[1].coordMat(surfQuadIndex,0);
        refCoord[1] = gaussPoints1d[quadIndex];
        // Evaluate all basis functions at point
        nodalBasis2d->evalBasis(refCoord,basisAtPoint);
        // Copy to recoveryInterpMatMu
        for (int basisIndex = 0; basisIndex < nlocal2d; basisIndex++)
          recoveryInterpMatMu[surfQuadIndex](quadIndex,basisIndex) = basisAtPoint[basisIndex];
      }
    }

    // Create a 2D region for (v,mu) looping for v
    int lower[2];
    int upper[2];
    lower[0] = localRgn.getLower(3);
    lower[1] = localRgn.getLower(4);
    upper[0] = localRgn.getUpper(3)+1;
    upper[1] = localRgn.getUpper(4);
    Lucee::Region<2, int> localRgn2d(lower,upper);
    // Loop over region
    //Lucee::RowMajorSequencer<2> seq(localRgn2d);
    Lucee::RowMajorSequencer<2> seqVPara(localRgn2d);
    Lucee::RowMajorIndexer<2> seqVParaIdxr(localRgn2d);
    int idx2d[2];
    int idx5dLower[5];
    int idx5dUpper[5];
    double xc[5];

    int totalSize = (localRgn2d.getUpper(0)-localRgn2d.getLower(0))*(localRgn2d.getUpper(1)-localRgn2d.getLower(1));
    recoveryMatricesVPara.resize(totalSize);

    for (int i = 0; i < 3; i++)
    {
      idx5dLower[i] = localRgn.getLower(i);
      idx5dUpper[i] = localRgn.getLower(i);
    }
    
    // Compute and store matrices for every cell
    while (seqVPara.step())
    {
      seqVPara.fillWithIndex(idx2d);
      int cellIndex = seqVParaIdxr.getIndex(idx2d);
      //std::cout << "Index " << cellIndex << " coords = [" << idx2d[0] << "," << idx2d[1] << "]" << std::endl;
      // At this interface (lower surface), compute recovery quadrature coordinates in v_parallel
      Eigen::VectorXd nodeList(2*numGaussPoints1d);
      idx5dLower[3] = idx2d[0]-1;
      idx5dLower[4] = idx2d[1];
      grid.setIndex(idx5dLower);
      grid.getCentroid(xc);
      for (int i = 0; i < numGaussPoints1d; i++)
        nodeList(i) = gaussPoints1d[i]*0.5*grid.getVolume()/grid.getSurfArea(3) - 0.5*grid.getVolume()/grid.getSurfArea(3);
      
      idx5dUpper[3] = idx2d[0];
      idx5dUpper[4] = idx2d[1];
      grid.setIndex(idx5dUpper);
      grid.getCentroid(xc);
      for (int i = numGaussPoints1d; i < 2*numGaussPoints1d; i++)
        nodeList(i) = gaussPoints1d[i-numGaussPoints1d]*0.5*grid.getVolume()/grid.getSurfArea(3) + 0.5*grid.getVolume()/grid.getSurfArea(3);
      //std::cout << "nodeList = " << nodeList << std::endl;

      Eigen::MatrixXd recoveryMatrix(2,2*numGaussPoints1d);
      for (int i = 0; i < recoveryMatrix.cols(); i++)
      {
        recoveryMatrix(0,i) = 1.0;
        recoveryMatrix(1,i) = 0.0;
        // Compute Lagrange interpolation polynomial evaluated at interface (x=0)
        for (int j = 0; j < recoveryMatrix.cols(); j++)
        {
          if (i != j)
          {
            recoveryMatrix(0,i) *= -nodeList(j)/(nodeList(i)-nodeList(j));
            recoveryMatrix(1,i) += 1/(-nodeList(j));
          }
        }
        recoveryMatrix(1,i) *= recoveryMatrix(0,i);
      }

      recoveryMatricesVPara[cellIndex] = recoveryMatrix;
    }

    // Do the same as above for the mu direction
    lower[0] = localRgn.getLower(3);
    lower[1] = localRgn.getLower(4);
    upper[0] = localRgn.getUpper(3);
    upper[1] = localRgn.getUpper(4)+1;
    Lucee::Region<2, int> localRgn2dMu(lower,upper);
    // Loop over region
    Lucee::RowMajorSequencer<2> seqMu(localRgn2dMu);
    Lucee::RowMajorIndexer<2> seqMuIdxr(localRgn2dMu);
    totalSize = (localRgn2dMu.getUpper(0)-localRgn2dMu.getLower(0))*(localRgn2dMu.getUpper(1)-localRgn2dMu.getLower(1));
    recoveryMatricesMu.resize(totalSize);
    // Compute and store matrices for every cell
    while (seqMu.step())
    {
      seqMu.fillWithIndex(idx2d);
      int cellIndex = seqMuIdxr.getIndex(idx2d);
      //std::cout << "Index " << cellIndex << " coords = [" << idx2d[0] << "," << idx2d[1] << "]" << std::endl;
      // At this interface (lower surface), compute recovery quadrature coordinates in mu
      Eigen::VectorXd nodeList(2*numGaussPoints1d);
      idx5dLower[3] = idx2d[0];
      idx5dLower[4] = idx2d[1]-1;
      grid.setIndex(idx5dLower);
      grid.getCentroid(xc);
      for (int i = 0; i < numGaussPoints1d; i++)
        nodeList(i) = gaussPoints1d[i]*0.5*grid.getVolume()/grid.getSurfArea(4) - 0.5*grid.getVolume()/grid.getSurfArea(4);
      
      idx5dUpper[3] = idx2d[0];
      idx5dUpper[4] = idx2d[1];
      grid.setIndex(idx5dUpper);
      grid.getCentroid(xc);
      for (int i = numGaussPoints1d; i < 2*numGaussPoints1d; i++)
        nodeList(i) = gaussPoints1d[i-numGaussPoints1d]*0.5*grid.getVolume()/grid.getSurfArea(4) + 0.5*grid.getVolume()/grid.getSurfArea(4);
      //std::cout << "nodeList = " << nodeList << std::endl;

      Eigen::MatrixXd recoveryMatrix(2,2*numGaussPoints1d);
      for (int i = 0; i < recoveryMatrix.cols(); i++)
      {
        recoveryMatrix(0,i) = 1.0;
        recoveryMatrix(1,i) = 0.0;
        // Compute Lagrange interpolation polynomial evaluated at interface (x=0)
        for (int j = 0; j < recoveryMatrix.cols(); j++)
        {
          if (i != j)
          {
            recoveryMatrix(0,i) *= -nodeList(j)/(nodeList(i)-nodeList(j));
            recoveryMatrix(1,i) += 1/(-nodeList(j));
          }
        }
        recoveryMatrix(1,i) *= recoveryMatrix(0,i);
      }

      recoveryMatricesMu[cellIndex] = recoveryMatrix;
    }
  }

  Lucee::UpdaterStatus
  LenardBernsteinDiffAlternate3D2VUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& fIn = this->getInp<Lucee::Field<5, double> >(0);
    // Temperature in joules
    const Lucee::Field<3, double>& temperatureIn = this->getInp<Lucee::Field<3, double> >(1);
    // Parallel Temperature in joules
    const Lucee::Field<3, double>& paraTemperatureIn = this->getInp<Lucee::Field<3, double> >(2);
    // Perpendicular Temperature in joules
    const Lucee::Field<3, double>& perpTemperatureIn = this->getInp<Lucee::Field<3, double> >(3);
    // Magnetic field
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(4);
    // Dimensionally correct number density from weighted moment calculation
    const Lucee::Field<3, double>& numDensityIn = this->getInp<Lucee::Field<3, double> >(5);
    // Output distribution function
    Lucee::Field<5, double>& fOut = this->getOut<Lucee::Field<5, double> >(0);

    double dt = t-this->getCurrTime();

    // Input fields
    Lucee::ConstFieldPtr<double> fInPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInLowerPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInUpperPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureInPtr = temperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> paraTemperatureInPtr = paraTemperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> perpTemperatureInPtr = perpTemperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> numDensityInPtr = numDensityIn.createConstPtr();
    // Writeable fields
    Lucee::FieldPtr<double> fOutPtr = fOut.createPtr();
    Lucee::FieldPtr<double> fOutLowerPtr = fOut.createPtr();
    Lucee::FieldPtr<double> fOutUpperPtr = fOut.createPtr();

    // check time-step
    double cflm = 1.1*cfl;
    double cfla = 0.0;
    
    fOut = 0.0;

    // local region to index
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<5> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[5];
    double cellCentroid[5];
    seq.fillWithIndex(idx);
    nodalBasis5d->setIndex(idx);
    int nlocal5d = nodalBasis5d->getNumNodes(); 
    nodalBasis3d->setIndex(idx[0], idx[1], idx[2]);
    int nlocal3d = nodalBasis3d->getNumNodes(); 
    nodalBasis2d->setIndex(idx[0], idx[1]);
    int nlocal2d = nodalBasis2d->getNumNodes(); 

    // Get alpha. Need to scale by n/T^(3/2) in this class
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    double alpha = resultVector[0];

    // Should be size 4 for linear elements
    Eigen::VectorXd fReduced(nodalStencil.size());
    Eigen::VectorXd fReducedLower(nodalStencil.size());
    Eigen::VectorXd fReducedUpper(nodalStencil.size());

    // Boundaries for surface integration
    /*int sliceLowerVPara = localRgn.getLower(3);
    int sliceUpperVPara = localRgn.getUpper(3)+1;
    // Make sure zero-flux BC's are enforced
    if (sliceLowerVPara == globalRgn.getLower(3))
      sliceLowerVPara = sliceLowerVPara + 1;
    if (sliceUpperVPara == globalRgn.getUpper(3)+1)
      sliceUpperVPara = sliceUpperVPara - 1;*/

    // Create a 2D region for looping and mapping surfaces in vPara (going from lower interior + 1 ghost cell)
    int lower[2];
    int upper[2];
    lower[0] = localRgn.getLower(3);
    lower[1] = localRgn.getLower(4);
    upper[0] = localRgn.getUpper(3)+1;
    upper[1] = localRgn.getUpper(4);
    Lucee::Region<2, int> localRgn2dVPara(lower,upper);
    Lucee::RowMajorSequencer<2> seqVPara(localRgn2dVPara);
    Lucee::RowMajorIndexer<2> seqVParaIdxr(localRgn2dVPara);
    // Create a 2D region for looping and mapping surfaces in mu (going from lower interior + 1 ghost cell)
    lower[0] = localRgn.getLower(3);
    lower[1] = localRgn.getLower(4);
    upper[0] = localRgn.getUpper(3);
    upper[1] = localRgn.getUpper(4)+1;
    Lucee::Region<2, int> localRgn2dMu(lower,upper);
    Lucee::RowMajorSequencer<2> seqMu(localRgn2dMu);
    Lucee::RowMajorIndexer<2> seqMuIdxr(localRgn2dMu);
    // Create a 2D region for looping over cells in mu (interior only)
    lower[0] = localRgn.getLower(3);
    lower[1] = localRgn.getLower(4);
    upper[0] = localRgn.getUpper(3);
    upper[1] = localRgn.getUpper(4);
    Lucee::Region<2, int> localRgn2dMuVol(lower,upper);
    Lucee::RowMajorSequencer<2> seqMuVol(localRgn2dMuVol);

    int idx2d[2];
    int idxLower[5];
    int idxUpper[5];

    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      idx[0] = ix;
      idxUpper[0] = ix;
      idxLower[0] = ix;
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        idx[1] = iy;
        idxUpper[1] = iy;
        idxLower[1] = iy;
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[2] = iz;
          idxUpper[2] = iz;
          idxLower[2] = iz;
          bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);
          temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
          paraTemperatureIn.setPtr(paraTemperatureInPtr, idx[0], idx[1], idx[2]);
          perpTemperatureIn.setPtr(perpTemperatureInPtr, idx[0], idx[1], idx[2]);
          numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);

          // At this location, loop over each configuration space grid node
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // V-parallel update
            seqVPara.reset();
            while (seqVPara.step())
            {
              seqVPara.fillWithIndex(idx2d);
              int cellIndex = seqVParaIdxr.getIndex(idx2d);
              // Check for zero-flux BC's
              if (idx2d[0] == globalRgn.getLower(3) || idx2d[0] == globalRgn.getUpper(3))
                continue;
              // If not a zero-flux BC cell, then proceed to compute recovery update on each side of surface
              idxLower[3] = idx2d[0]-1;
              idxUpper[3] = idx2d[0];
              idxLower[4] = idx2d[1];
              idxUpper[4] = idx2d[1];

              fIn.setPtr(fInLowerPtr, idxLower);
              fIn.setPtr(fInUpperPtr, idxUpper);
              
              // Fill out fReduced at this location
              for (int i = 0; i < fReduced.size(); i++)
              {
                fReducedLower(i) = fInLowerPtr[configNode + nodalStencil[i]];
                fReducedUpper(i) = fInUpperPtr[configNode + nodalStencil[i]];
              }
              
              // Matrix storing recovery polynomial and its derivative at each surface quadrature point
              // recoveryInterMap.size() should equal number of surface quadrature points
              Eigen::MatrixXd recoveryEval = Eigen::MatrixXd::Zero(recoveryInterpMatVPara.size(), 2);
              // Loop over each surface quadrature point
              for (int surfQuadIndex = 0; surfQuadIndex < recoveryInterpMatVPara.size(); surfQuadIndex++)
              {
                // Compute f at recovery quadrature points
                Eigen::VectorXd fLowerAtRecoveryQuad = recoveryInterpMatVPara[surfQuadIndex]*fReducedLower;
                Eigen::VectorXd fUpperAtRecoveryQuad = recoveryInterpMatVPara[surfQuadIndex]*fReducedUpper;

                // Join together fLowerAtRecoveryQuad and fUpperAtRecoveryQuad
                Eigen::VectorXd fRecoveryQuad(2*fLowerAtRecoveryQuad.size());
                fRecoveryQuad << fLowerAtRecoveryQuad,fUpperAtRecoveryQuad;
                
                // Compute recovery polynomial and its derivative at this surface quadrature point
                recoveryEval.row(surfQuadIndex) = recoveryMatricesVPara[cellIndex]*fRecoveryQuad;
                //std::cout << "recoverySolution" << std::endl << recoveryMatricesVPara[cellIndex]*fRecoveryQuad << std::endl;
                //std::cout << "recoveryEval" << std::endl << recoveryEval << std::endl;
              }

              // With recovery polynomial and its derivative known at every surface quadrature point, ready
              // to perform updates to distribution function in lower and upper cells
              grid.setIndex(idxUpper);
              double term1Scale = grid.getSurfArea(3)/grid.getVolume();
              double term2Scale = grid.getSurfArea(3)/grid.getVolume()*grid.getSurfArea(3)/grid.getVolume();
              Eigen::VectorXd upperResultVector = term1Scale*surfLowerUpdateMatrix[0]*recoveryEval.col(1) - 
                term2Scale*surfLowerUpdateMatrixDeriv[0]*recoveryEval.col(0);
              // Keep track of max CFL number
              cfla = std::max( cfla, std::abs(4.0*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))*
                paraTemperatureInPtr[configNode]/speciesMass*dt/(grid.getVolume()/grid.getSurfArea(3)*grid.getVolume()/grid.getSurfArea(3))) );

              grid.setIndex(idxLower);
              term1Scale = grid.getSurfArea(3)/grid.getVolume();
              term2Scale = grid.getSurfArea(3)/grid.getVolume()*grid.getSurfArea(3)/grid.getVolume();
              Eigen::VectorXd lowerResultVector = term1Scale*surfUpperUpdateMatrix[0]*recoveryEval.col(1) - 
                term2Scale*surfUpperUpdateMatrixDeriv[0]*recoveryEval.col(0);
              // Keep track of max CFL number
              cfla = std::max( cfla, std::abs(4.0*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))*
                paraTemperatureInPtr[configNode]/speciesMass*dt/(grid.getVolume()/grid.getSurfArea(3)*grid.getVolume()/grid.getSurfArea(3))) );

              // Accumulate to output
              fOut.setPtr(fOutLowerPtr, idxLower);
              fOut.setPtr(fOutUpperPtr, idxUpper);
              for (int i = 0; i < fReduced.size(); i++)
              {
                fOutLowerPtr[configNode + nodalStencil[i]] += paraTemperatureInPtr[configNode]/speciesMass*lowerResultVector(i);
                fOutUpperPtr[configNode + nodalStencil[i]] -= paraTemperatureInPtr[configNode]/speciesMass*upperResultVector(i);
              }
            }

            // Mu update (Volumes)
            seqMuVol.reset();
            while (seqMuVol.step())
            {
              seqMuVol.fillWithIndex(idx2d);
              
              idx[3] = idx2d[0];
              idx[4] = idx2d[1];

              fIn.setPtr(fInPtr, idx);

              // Fill out fReduced at this location
              for (int i = 0; i < fReduced.size(); i++)
                fReduced(i) = fInPtr[configNode + nodalStencil[i]];
              
              grid.setIndex(idx);
              // Compute volume update (with 1/dmu scale factor for grad stiffness matrix)
              Eigen::VectorXd resultVector = 2*perpTemperatureInPtr[configNode]/bFieldInPtr[configNode]*
                grid.getSurfArea(4)/grid.getVolume()*volumeUpdateMatrix[1]*fReduced;

              // CFL number check
              grid.setIndex(idx);
              grid.getCentroid(cellCentroid);

              double muTherm = perpTemperatureInPtr[configNode]/bFieldInPtr[configNode];
              
              for (int quadIndex = 0; quadIndex < gaussPoints1d.size(); quadIndex++)
              {
                double muCoord = cellCentroid[4] + 0.5*grid.getVolume()/grid.getSurfArea(4)*gaussPoints1d[quadIndex];
                cfla = std::max(cfla, 8.0*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))
                *muTherm*muCoord*dt/(grid.getVolume()/grid.getSurfArea(4)*grid.getVolume()/grid.getSurfArea(4)));
              }

              // Accumulate to output
              fOut.setPtr(fOutPtr, idx);
              for (int i = 0; i < fReduced.size(); i++)
                fOutPtr[configNode + nodalStencil[i]] += resultVector(i);
            }

            // Mu update (Surfaces)
            seqMu.reset();
            while (seqMu.step())
            {
              seqMu.fillWithIndex(idx2d);
              int cellIndex = seqMuIdxr.getIndex(idx2d);
              // Check for zero-flux BC's
              if (idx2d[1] == globalRgn.getLower(4) || idx2d[1] == globalRgn.getUpper(4))
                continue;
              // If not a zero-flux BC cell, then proceed to compute recovery update on each side of surface
              idxLower[3] = idx2d[0];
              idxUpper[3] = idx2d[0];
              idxLower[4] = idx2d[1]-1;
              idxUpper[4] = idx2d[1];

              fIn.setPtr(fInLowerPtr, idxLower);
              fIn.setPtr(fInUpperPtr, idxUpper);
              
              // Fill out fReduced at this location
              for (int i = 0; i < fReduced.size(); i++)
              {
                fReducedLower(i) = fInLowerPtr[configNode + nodalStencil[i]];
                fReducedUpper(i) = fInUpperPtr[configNode + nodalStencil[i]];
              }
              
              // Matrix storing recovery polynomial and its derivative at each surface quadrature point
              // recoveryInterMap.size() should equal number of surface quadrature points
              Eigen::MatrixXd recoveryEval = Eigen::MatrixXd::Zero(recoveryInterpMatMu.size(), 2);
              // Loop over each surface quadrature point
              for (int surfQuadIndex = 0; surfQuadIndex < recoveryInterpMatMu.size(); surfQuadIndex++)
              {
                // Compute f at recovery quadrature points
                Eigen::VectorXd fLowerAtRecoveryQuad = recoveryInterpMatMu[surfQuadIndex]*fReducedLower;
                Eigen::VectorXd fUpperAtRecoveryQuad = recoveryInterpMatMu[surfQuadIndex]*fReducedUpper;

                // Join together fLowerAtRecoveryQuad and fUpperAtRecoveryQuad
                Eigen::VectorXd fRecoveryQuad(2*fLowerAtRecoveryQuad.size());
                fRecoveryQuad << fLowerAtRecoveryQuad,fUpperAtRecoveryQuad;
                
                // Compute recovery polynomial and its derivative at this surface quadrature point
                recoveryEval.row(surfQuadIndex) = recoveryMatricesMu[cellIndex]*fRecoveryQuad;
                //std::cout << "recoverySolution" << std::endl << recoveryMatricesVPara[cellIndex]*fRecoveryQuad << std::endl;
                //std::cout << "recoveryEval" << std::endl << recoveryEval << std::endl;
              }
              
              // With recovery polynomial and its derivative known at every surface quadrature point, ready
              // to perform updates to distribution function in lower and upper cells
              grid.setIndex(idxUpper);
              grid.getCentroid(cellCentroid);
              double muCoord = cellCentroid[4] - 0.5*grid.getVolume()/grid.getSurfArea(4);
              
              double term1Scale = muCoord*grid.getSurfArea(4)/grid.getVolume();
              double term2Scale = muCoord*grid.getSurfArea(4)/grid.getVolume()*grid.getSurfArea(4)/grid.getVolume();
              Eigen::VectorXd upperResultVector = term1Scale*surfLowerUpdateMatrix[1]*recoveryEval.col(1) - 
                term2Scale*surfLowerUpdateMatrixDeriv[1]*recoveryEval.col(0);

              grid.setIndex(idxLower);
              term1Scale = muCoord*grid.getSurfArea(4)/grid.getVolume();
              term2Scale = muCoord*grid.getSurfArea(4)/grid.getVolume()*grid.getSurfArea(4)/grid.getVolume();
              Eigen::VectorXd lowerResultVector = term1Scale*surfUpperUpdateMatrix[1]*recoveryEval.col(1) - 
                term2Scale*surfUpperUpdateMatrixDeriv[1]*recoveryEval.col(0);

              // Accumulate to output
              fOut.setPtr(fOutLowerPtr, idxLower);
              fOut.setPtr(fOutUpperPtr, idxUpper);
              for (int i = 0; i < fReduced.size(); i++)
              {
                fOutLowerPtr[configNode + nodalStencil[i]] += 2*perpTemperatureInPtr[configNode]/bFieldInPtr[configNode]*lowerResultVector(i);
                fOutUpperPtr[configNode + nodalStencil[i]] -= 2*perpTemperatureInPtr[configNode]/bFieldInPtr[configNode]*upperResultVector(i);
              }
            }
          }
        }
      }
    }

    seq.reset();
    // Final sweep, update solution with forward Euler step
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fOut.setPtr(fOutPtr, idx);
      // Compute alpha scale factor n/Te^(3/2)
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);

      if (onlyIncrement == false)
      {
        fIn.setPtr(fInPtr, idx);
        for (int configNode = 0; configNode < nlocal3d; configNode++)
        {
          for (int stencilIndex = 0; stencilIndex < nodalStencil.size(); stencilIndex++)
          {
            fOutPtr[configNode + nodalStencil[stencilIndex]] = fInPtr[configNode + nodalStencil[stencilIndex]]
              + dt*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))
                  *fOutPtr[configNode + nodalStencil[stencilIndex]];
          }
        }
      }
      else
      {
        for (int configNode = 0; configNode < nlocal3d; configNode++)
        {
          for (int stencilIndex = 0; stencilIndex < nodalStencil.size(); stencilIndex++)
          {
            fOutPtr[configNode + nodalStencil[stencilIndex]] = alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]
                *sqrt(temperatureInPtr[configNode]))
                *fOutPtr[configNode + nodalStencil[stencilIndex]];
          }
        }
      }
    }
    
    if (cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);
    else
      return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDiffAlternate3D2VUpdater::declareTypes()
  {
    // Input: Distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Parallel Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Perpendicular Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Magnetic field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: dimensionally correct number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // returns one output (fNew)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  LenardBernsteinDiffAlternate3D2VUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  LenardBernsteinDiffAlternate3D2VUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("LenardBernsteinDiffAlternate3D2VUpdater::evaluateFunction: ");
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
        throw Lucee::Except("LenardBernsteinDiffAlternate3D2VUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

  bool
  LenardBernsteinDiffAlternate3D2VUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

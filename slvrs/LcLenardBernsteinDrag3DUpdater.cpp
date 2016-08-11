/**
 * @file	LcLenardBernsteinDrag3DUpdater.cpp
 *
 * @brief	Updater to compute the drag term in the L-B collision operator.
 * df/dt = alpha*d[(v-u)f]/dv where u is a function of x (drift velocity)
 * Used for 1D2V SOL problem.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLenardBernsteinDrag3DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *LenardBernsteinDrag3DUpdater::id = "LenardBernsteinDrag3DUpdater";

  LenardBernsteinDrag3DUpdater::LenardBernsteinDrag3DUpdater()
    : UpdaterIfc()
  {
  }

  void 
  LenardBernsteinDrag3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("LenardBernsteinDrag3DUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except("LenardBernsteinDrag3DUpdater::readInput: Must specify element to use using 'basis1d'");
    
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
      throw Lucee::Except("LenardBernsteinDrag3DUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void 
  LenardBernsteinDrag3DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();
    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<3> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[3];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    nodalBasis1d->setIndex(idx[0]);
    
    int nlocal = nodalBasis->getNumNodes();
    int nlocal1d = nodalBasis1d->getNumNodes();

    int numVolQuadNodes = nodalBasis->getNumGaussNodes();
    int numSurfQuadNodes = nodalBasis->getNumSurfGaussNodes();

    int nVolQuad1d = nodalBasis1d->getNumGaussNodes();
    int numVolQuadNodes1d = nodalBasis1d->getNumGaussNodes();

    // Get mass and grad-stiffness matrices and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);
    Lucee::Matrix<double> gradStiffMatrixLucee(nlocal, nlocal);

    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    // Element 0 is dir 1, Element 1 is dir 2
    std::vector<Eigen::MatrixXd> gradStiffMatrix(2, Eigen::MatrixXd(nlocal, nlocal));
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    
    nodalBasis->getGradStiffnessMatrix(1, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix[0]);
    nodalBasis->getGradStiffnessMatrix(2, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix[1]);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpVolMatrixLucee(numVolQuadNodes, nlocal);
    Lucee::Matrix<double> interpSurfMatrixLowerLucee(numSurfQuadNodes, nlocal);
    Lucee::Matrix<double> interpSurfMatrixUpperLucee(numSurfQuadNodes, nlocal);
    Lucee::Matrix<double> interpVolMatrix1dLucee(numVolQuadNodes1d, nlocal1d);

    Lucee::Matrix<double> gaussVolOrdinatesLucee(numVolQuadNodes, 3);
    Lucee::Matrix<double> gaussVolOrdinates1dLucee(numVolQuadNodes1d, 3);
    Lucee::Matrix<double> gaussSurfOrdinatesLucee(numSurfQuadNodes, 3);

    gaussVolWeights = std::vector<double>(numVolQuadNodes);
    std::vector<double> gaussVolWeights1d(numVolQuadNodes1d);
    gaussVolOrdinates = Eigen::MatrixXd(numVolQuadNodes, 3);
    // Get the interpolation matrix for the volume quadrature points
    interpVolMatrix = Eigen::MatrixXd(numVolQuadNodes, nlocal);
    nodalBasis->getGaussQuadData(interpVolMatrixLucee, gaussVolOrdinatesLucee, gaussVolWeights);
    copyLuceeToEigen(interpVolMatrixLucee, interpVolMatrix);
    copyLuceeToEigen(gaussVolOrdinatesLucee, gaussVolOrdinates);
    // Get the 1d interpolation matrix
    interpVolMatrix1d = Eigen::MatrixXd(numVolQuadNodes1d, nlocal1d);
    nodalBasis1d->getGaussQuadData(interpVolMatrix1dLucee, gaussVolOrdinates1dLucee, gaussVolWeights1d);
    copyLuceeToEigen(interpVolMatrix1dLucee, interpVolMatrix1d);

    interpSurfMatrixLower.resize(2, Eigen::MatrixXd(numSurfQuadNodes, nlocal));
    interpSurfMatrixUpper.resize(2, Eigen::MatrixXd(numSurfQuadNodes, nlocal));
    gaussSurfOrdinates.resize(2, Eigen::MatrixXd(numSurfQuadNodes, 3));
    gaussSurfWeights.resize(2, std::vector<double>(numSurfQuadNodes));
    // Get the interpolation matrix for the upper surface quadrature points.
    // Quadrature location and weights will be overwritten but it doesn't matter for these purposes
    nodalBasis->getSurfUpperGaussQuadData(1, interpSurfMatrixUpperLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights[0]);
    copyLuceeToEigen(interpSurfMatrixUpperLucee, interpSurfMatrixUpper[0]);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis->getSurfLowerGaussQuadData(1, interpSurfMatrixLowerLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights[0]);
    copyLuceeToEigen(interpSurfMatrixLowerLucee, interpSurfMatrixLower[0]);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates[0]);

    // Do the same thing as above in the mu direction
    nodalBasis->getSurfUpperGaussQuadData(2, interpSurfMatrixUpperLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights[1]);
    copyLuceeToEigen(interpSurfMatrixUpperLucee, interpSurfMatrixUpper[1]);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis->getSurfLowerGaussQuadData(2, interpSurfMatrixLowerLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights[1]);
    copyLuceeToEigen(interpSurfMatrixLowerLucee, interpSurfMatrixLower[1]);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates[1]);

    // Compute and store inverse of mass matrix
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();
    // Compute derivative of basis function evaluated at volume quadrature points
    basisDerivAtVolQuad.resize(2);
    basisDerivAtVolQuad[0] = interpVolMatrix*massMatrixInv*gradStiffMatrix[0].transpose();
    basisDerivAtVolQuad[1] = interpVolMatrix*massMatrixInv*gradStiffMatrix[1].transpose();
    // Take transpose so same basis function evaluated at diff quad points in the same row
    // Need to use transposeInPlace because otherwise there will be a bug!
    basisDerivAtVolQuad[0].transposeInPlace();
    basisDerivAtVolQuad[1].transposeInPlace();

    // Pre-multiply by inverse mass matrix
    basisDerivAtVolQuad[0] = massMatrixInv*basisDerivAtVolQuad[0];
    basisDerivAtVolQuad[1] = massMatrixInv*basisDerivAtVolQuad[1];

    // Precompute two matrices for surface integrals
    surfIntegralMatrixLower.resize(2);
    surfIntegralMatrixUpper.resize(2);
    // Pre-multiply by inverse mass matrix
    surfIntegralMatrixLower[0] = massMatrixInv*interpSurfMatrixLower[0].transpose();
    surfIntegralMatrixUpper[0] = massMatrixInv*interpSurfMatrixUpper[0].transpose();
    surfIntegralMatrixLower[1] = massMatrixInv*interpSurfMatrixLower[1].transpose();
    surfIntegralMatrixUpper[1] = massMatrixInv*interpSurfMatrixUpper[1].transpose();
  }

  Lucee::UpdaterStatus 
  LenardBernsteinDrag3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Distribution function
    const Lucee::Field<3, double>& q = this->getInp<Lucee::Field<3, double> >(0);
    // Drift velocity u(x).
    const Lucee::Field<1, double>& u = this->getInp<Lucee::Field<1, double> >(1);
    // Output distribution function
    Lucee::Field<3, double>& qNew = this->getOut<Lucee::Field<3, double> >(0);

    int nlocal = nodalBasis->getNumNodes();
    int nlocal1d = nodalBasis1d->getNumNodes();
    int nVolQuad1d = nodalBasis1d->getNumGaussNodes();

    double dt = t-this->getCurrTime();

    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();

    double cfla = 0.0; // maximum CFL number

    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::ConstFieldPtr<double> uPtr = u.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrLeft = qNew.createPtr();

    qNew = 0.0; // use qNew to store increment initially
    int idx[3];
    double cellCentroid[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);
    
    // Get alpha
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    alpha = resultVector[0];

    // Volume integral contribution
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(qPtr, idx);
      qNew.setPtr(qNewPtr, idx);
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);
      
      u.setPtr(uPtr, idx[0]);
      Eigen::VectorXd uVals(nlocal1d);
        
      // Copy values of u into uVals vector
      for (int nodeIndex = 0; nodeIndex < nlocal1d; nodeIndex++)
        uVals(nodeIndex) = uPtr[nodeIndex];
      // Interpolate u to 1D quadrature points 
      Eigen::VectorXd uAtQuad = interpVolMatrix1d*uVals;
        
      // Figure out what f is at each quadrature point
      Eigen::VectorXd fVals(nlocal);
      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        fVals(nodeIndex) = qPtr[nodeIndex];
      
      Eigen::VectorXd fVolQuad = interpVolMatrix*fVals;

      Eigen::VectorXd paraQuad(fVolQuad.size());
      Eigen::VectorXd muQuad(fVolQuad.size());
      // Compute integrands for gaussian quadrature excluding basis function derivatives
      for (int quadIndex = 0; quadIndex < fVolQuad.size(); quadIndex++)
      {
        double paraCoord = cellCentroid[1] + 0.5*grid.getDx(1)*gaussVolOrdinates(quadIndex,1);
        double muCoord = cellCentroid[2] + 0.5*grid.getDx(2)*gaussVolOrdinates(quadIndex,2);
        paraQuad(quadIndex) = gaussVolWeights[quadIndex]*fVolQuad(quadIndex)*
          (paraCoord - uAtQuad(quadIndex % nVolQuad1d));
        muQuad(quadIndex) = gaussVolWeights[quadIndex]*fVolQuad(quadIndex)*2*muCoord;
      }

      // Evaluate integral using gaussian quadrature (represented as matrix-vector multiply)
      Eigen::VectorXd volIntegralResult = (basisDerivAtVolQuad[0]*paraQuad +
        basisDerivAtVolQuad[1]*muQuad);
      //Eigen::VectorXd volIntegralResult = basisDerivAtVolQuad[1]*muQuad;
          
      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        qNewPtr[nodeIndex] -= volIntegralResult(nodeIndex);
    }

    // Contributions from surface integrals
    // Create sequencer to loop over *each* 2D slice in '1' direction
    Lucee::RowMajorSequencer<3> seqLowerDim(localRgn.deflate(1));

    int idxr[3], idxl[3];
    
    int ivLower = localRgn.getLower(1);
    // Need one edge outside domain interior
    int ivUpper = localRgn.getUpper(1)+1;

    // Zero flux BC's
    if (ivLower == globalRgn.getLower(1))
      ivLower = globalRgn.getLower(1)+1;
    if (ivUpper == globalRgn.getUpper(1)+1)
      ivUpper = globalRgn.getUpper(1);

    // Loop over each 1D slice in vParallel
    while (seqLowerDim.step())
    {
      seqLowerDim.fillWithIndex(idxr);
      seqLowerDim.fillWithIndex(idxl);

      // Compute u at quadrature points (again)
      u.setPtr(uPtr, idxr[0]);
      Eigen::VectorXd uAtNodes(nlocal1d);
      for(int nodeIndex = 0; nodeIndex < nlocal1d; nodeIndex++)
        uAtNodes(nodeIndex) = uPtr[nodeIndex];
      Eigen::VectorXd uAtQuad = interpVolMatrix1d*uAtNodes;

      // Loop over each edge in vParallel
      for (int i = ivLower; i < ivUpper; i++)
      {
        idxr[1] = i; // cell right of edge
        idxl[1] = i-1; // cell left of edge
        
        q.setPtr(qPtr, idxr);
        q.setPtr(qPtrl, idxl);
        
        grid.setIndex(idxl);
        double dxL = grid.getDx(1);
        grid.setIndex(idxr);
        double dxR = grid.getDx(1);
        // Need to know center of right cell to figure out
        // the global velocity coordinate
        grid.getCentroid(cellCentroid);

        // Copy q data to Eigen vectors for matrix multiplications
        Eigen::VectorXd fLeft(nlocal);
        Eigen::VectorXd fRight(nlocal);
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          fLeft(nodeIndex) = qPtrl[nodeIndex];
          fRight(nodeIndex) = qPtr[nodeIndex];
        }

        // Evaluate fLeft and fRight at quadrature nodes on surface
        Eigen::VectorXd fLeftSurfEvals = interpSurfMatrixUpper[0]*fLeft;
        Eigen::VectorXd fRightSurfEvals = interpSurfMatrixLower[0]*fRight;

        // Loop over quadrature points on the edge
        Eigen::VectorXd surfIntegralFluxes(gaussSurfOrdinates[0].rows());

        double paraCoord = cellCentroid[1] + gaussSurfOrdinates[0](0,1)*0.5*grid.getDx(1);
        for (int quadIndex = 0; quadIndex < gaussSurfOrdinates[0].rows(); quadIndex++)
        {

          // Compute Lax flux at each surface quadrature point
          double numFlux;
          if (paraCoord - uAtQuad(quadIndex % nVolQuad1d) > 0)
            numFlux = (paraCoord - uAtQuad(quadIndex % nVolQuad1d))*fRightSurfEvals(quadIndex);
          else numFlux = (paraCoord - uAtQuad(quadIndex % nVolQuad1d))*fLeftSurfEvals(quadIndex);

          //numFlux = 0.5*((paraCoord - uAtQuad(quadIndex % nVolQuad1d))*fRightSurfEvals(quadIndex) +
          //  (paraCoord - uAtQuad(quadIndex % nVolQuad1d))*fLeftSurfEvals(quadIndex));
          // Store result of weight*(v-u)*f at this location in a vector
          surfIntegralFluxes(quadIndex) = gaussSurfWeights[0][quadIndex]*numFlux;

          // Keep track of max CFL number
          cfla = std::max(cfla, std::abs(alpha*(paraCoord-uAtQuad(quadIndex % nVolQuad1d))*dt/grid.getDx(1)));
          //cfla = std::max(cfla, std::abs(alpha*(paraCoord-uAtQuad(quadIndex % nVolQuad1d))*dt/(0.5*(dxL+dxR))));
          // Time-step was too large: return a suggestion with correct time-step
          if (cfla > cflm)
            return Lucee::UpdaterStatus(false, dt*cfl/cfla);
        }

        // Compute all surface integrals using a matrix multiply.
        // Each row in result is a basis function times flux projection
        // Then inverse of mass matrix is multiplied to find appropriate increments
        Eigen::VectorXd leftSurfIntegralResult = surfIntegralMatrixUpper[0]*surfIntegralFluxes;
        Eigen::VectorXd rightSurfIntegralResult = surfIntegralMatrixLower[0]*surfIntegralFluxes;

        // Update left cell connected to edge with flux on face
        qNew.setPtr(qNewPtrLeft, idxl);
        // Update right cell connected to edge with flux on face
        qNew.setPtr(qNewPtr, idxr);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          // Plus sign for left cell since outward normal is in +v direction
          qNewPtrLeft[nodeIndex] += leftSurfIntegralResult(nodeIndex);
          // Minus sign for right cell since outward normal is in -v direction
          qNewPtr[nodeIndex]  -= rightSurfIntegralResult(nodeIndex);
        }
      }
    }

    // Contributions from surface integrals
    // Create sequencer to loop over *each* 2D slice in '2' direction
    seqLowerDim = Lucee::RowMajorSequencer<3>(localRgn.deflate(2));

    int iMuLower = localRgn.getLower(2);
    // Need one edge outside domain interior
    int iMuUpper = localRgn.getUpper(2)+1;

    // Zero flux BC's
    if (iMuLower == globalRgn.getLower(2))
      iMuLower = globalRgn.getLower(2)+1;
    if (iMuUpper == globalRgn.getUpper(2)+1)
      iMuUpper = globalRgn.getUpper(2);
  
    // Loop over each 1D slice in mu
    while (seqLowerDim.step())
    {
      seqLowerDim.fillWithIndex(idxr);
      seqLowerDim.fillWithIndex(idxl);

      // Loop over each edge in mu
      for (int i = iMuLower; i < iMuUpper; i++)
      {
        idxr[2] = i; // cell right of edge
        idxl[2] = i-1; // cell left of edge
        
        q.setPtr(qPtr, idxr);
        q.setPtr(qPtrl, idxl);
        
        grid.setIndex(idxl);
        double dxL = grid.getDx(2);
        grid.setIndex(idxr);
        double dxR = grid.getDx(2);
        // Need to know center of right cell to figure out
        // the global velocity coordinate
        grid.getCentroid(cellCentroid);

        // Copy q data to Eigen vectors for matrix multiplications
        Eigen::VectorXd fLeft(nlocal);
        Eigen::VectorXd fRight(nlocal);
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          fLeft(nodeIndex) = qPtrl[nodeIndex];
          fRight(nodeIndex) = qPtr[nodeIndex];
        }

        // Evaluate fLeft and fRight at quadrature nodes on surface
        Eigen::VectorXd fLeftSurfEvals = interpSurfMatrixUpper[1]*fLeft;
        Eigen::VectorXd fRightSurfEvals = interpSurfMatrixLower[1]*fRight;

        // Loop over quadrature points on the edge
        Eigen::VectorXd surfIntegralFluxes(gaussSurfOrdinates[1].rows());

        double muCoord = cellCentroid[2] + gaussSurfOrdinates[1](0,2)*0.5*grid.getDx(2);
        for (int quadIndex = 0; quadIndex < gaussSurfOrdinates[1].rows(); quadIndex++)
        {
          // Compute Lax flux at each surface quadrature point
          double numFlux;
          if (muCoord > 0)
            numFlux = 2*muCoord*fRightSurfEvals(quadIndex);
          else numFlux = 2*muCoord*fLeftSurfEvals(quadIndex);

          //numFlux = 0.5*(2*muCoord*fRightSurfEvals(quadIndex)+2*muCoord*fLeftSurfEvals(quadIndex));
          // Store result of weight*mu*f at this location in a vector
          surfIntegralFluxes(quadIndex) = gaussSurfWeights[1][quadIndex]*numFlux;

          // Keep track of max CFL number
          cfla = std::max(cfla, std::abs(alpha*2*muCoord*dt/grid.getDx(2)));
          // Time-step was too large: return a suggestion with correct time-step
          if (cfla > cflm)
            return Lucee::UpdaterStatus(false, dt*cfl/cfla);
        }

        // Compute all surface integrals using a matrix multiply.
        // Each row in result is a basis function times flux projection
        // Then inverse of mass matrix is multiplied to find appropriate increments
        Eigen::VectorXd leftSurfIntegralResult = surfIntegralMatrixUpper[1]*surfIntegralFluxes;
        Eigen::VectorXd rightSurfIntegralResult = surfIntegralMatrixLower[1]*surfIntegralFluxes;

        // Update left cell connected to edge with flux on face
        qNew.setPtr(qNewPtrLeft, idxl);
        // Update right cell connected to edge with flux on face
        qNew.setPtr(qNewPtr, idxr);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          // Plus sign for left cell since outward normal is in +v direction
          qNewPtrLeft[nodeIndex] += leftSurfIntegralResult(nodeIndex);
          // Minus sign for right cell since outward normal is in -v direction
          qNewPtr[nodeIndex] -= rightSurfIntegralResult(nodeIndex);
        }
      }
    }
    
    seq.reset();// = Lucee::RowMajorSequencer<3>(localRgn);
    // Final sweep, update solution with forward Euler step

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      qNew.setPtr(qNewPtr, idx);
      
      if (onlyIncrement == false)
      {
        q.setPtr(qPtr, idx);
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          qNewPtr[nodeIndex] = qPtr[nodeIndex] + dt*alpha*qNewPtr[nodeIndex];
      }
      else
      {
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          qNewPtr[nodeIndex] = alpha*qNewPtr[nodeIndex];
      }
    }
    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDrag3DUpdater::declareTypes()
  {
    // takes two inputs (fOld, u_parallel) 
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output (fNew)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  LenardBernsteinDrag3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  LenardBernsteinDrag3DUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("HeatFluxAtEdgeUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'tPerpProfile' "
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
        throw Lucee::Except("HeatFluxAtEdgeUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}

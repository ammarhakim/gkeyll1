/**
 * @file	LcLenardBernsteinDragPara3D2VUpdater.cpp
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
#include <LcLenardBernsteinDragPara3D2VUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *LenardBernsteinDragPara3D2VUpdater::id = "LenardBernsteinDragPara3D2VUpdater";

  LenardBernsteinDragPara3D2VUpdater::LenardBernsteinDragPara3D2VUpdater()
    : UpdaterIfc()
  {
  }

  void 
  LenardBernsteinDragPara3D2VUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("LenardBernsteinDragPara3D2VUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("LenardBernsteinDragPara3D2VUpdater::readInput: Must specify element to use using 'basis3d'");
    
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("LenardBernsteinDragPara3D2VUpdater::readInput: Must specify element to use using 'basis2d'");
    
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
      throw Lucee::Except("LenardBernsteinDragPara3D2VUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void 
  LenardBernsteinDragPara3D2VUpdater::initialize()
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
    nodalBasis2d->setIndex(idx[0], idx[1]);
    nodalBasis5d->setIndex(idx);

    int nlocal5d = nodalBasis5d->getNumNodes();
    int nlocal2d = nodalBasis2d->getNumNodes();
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    int nSurfQuad2d = nodalBasis2d->getNumSurfGaussNodes();

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

    // Get mass and grad-stiffness matrices and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal2d, nlocal2d);
    Lucee::Matrix<double> gradStiffMatrixLucee(nlocal2d, nlocal2d);

    Eigen::MatrixXd massMatrix(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(massMatrixLucee);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    massMatrix *= grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));

    // Element 0 is vpar (dir 0 in 2d),  Element 1 is mu (dir 1 in 2d)
    std::vector<Eigen::MatrixXd> gradStiffMatrix(2, Eigen::MatrixXd(nlocal2d, nlocal2d));
    nodalBasis2d->getGradStiffnessMatrix(0, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix[0]);
    gradStiffMatrix[0] *= grid.getDx(4)/grid.getDx(1);

    nodalBasis2d->getGradStiffnessMatrix(1, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix[1]);
    gradStiffMatrix[1] *= grid.getDx(3)/grid.getDx(0);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpVolMatrixLucee(nVolQuad2d, nlocal2d);
    Lucee::Matrix<double> interpSurfMatrixLower2dLucee(nSurfQuad2d, nlocal2d);
    Lucee::Matrix<double> interpSurfMatrixUpper2dLucee(nSurfQuad2d, nlocal2d);

    Lucee::Matrix<double> gaussVolOrdinatesLucee(nVolQuad2d, 3);
    Lucee::Matrix<double> gaussSurfOrdinatesLucee(nSurfQuad2d, 3);
    gaussVolWeights2d = std::vector<double>(nVolQuad2d);
    gaussVolOrdinates2d = Eigen::MatrixXd(nVolQuad2d, 3);
    
    // Get the interpolation matrix for the volume quadrature points
    interpVolMatrix2d = Eigen::MatrixXd(nVolQuad2d, nlocal2d);
    nodalBasis2d->getGaussQuadData(interpVolMatrixLucee, gaussVolOrdinatesLucee, gaussVolWeights2d);

    // Scale gaussVolWeights2d for (v,mu) integration
    double scaleCorrection = grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));
    for (int i = 0; i < gaussVolWeights2d.size(); i++)
      gaussVolWeights2d[i] = scaleCorrection*gaussVolWeights2d[i];
    
    copyLuceeToEigen(interpVolMatrixLucee, interpVolMatrix2d);
    copyLuceeToEigen(gaussVolOrdinatesLucee, gaussVolOrdinates2d);

    interpSurfMatrixLower2d.resize(2, Eigen::MatrixXd(nSurfQuad2d, nlocal2d));
    interpSurfMatrixUpper2d.resize(2, Eigen::MatrixXd(nSurfQuad2d, nlocal2d));
    gaussSurfOrdinates2d.resize(2, Eigen::MatrixXd(nSurfQuad2d, 3));
    gaussSurfWeights2d.resize(2, std::vector<double>(nSurfQuad2d));

    // Get the interpolation matrix for the upper surface quadrature points.
    // Quadrature location and weights will be overwritten but it doesn't matter for these purposes
    nodalBasis2d->getSurfUpperGaussQuadData(0, interpSurfMatrixUpper2dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights2d[0]);
    copyLuceeToEigen(interpSurfMatrixUpper2dLucee, interpSurfMatrixUpper2d[0]);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis2d->getSurfLowerGaussQuadData(0, interpSurfMatrixLower2dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights2d[0]);
    copyLuceeToEigen(interpSurfMatrixLower2dLucee, interpSurfMatrixLower2d[0]);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates2d[0]);
    
    // Scale gaussSurfWeights2d[0] for (v,mu) integration
    for (int i = 0; i < gaussSurfWeights2d[0].size(); i++)
      gaussSurfWeights2d[0][i] = gaussSurfWeights2d[0][i]*grid.getDx(4)/grid.getDx(1);

    // Do the same thing as above in the mu direction
    nodalBasis2d->getSurfUpperGaussQuadData(1, interpSurfMatrixUpper2dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights2d[1]);
    copyLuceeToEigen(interpSurfMatrixUpper2dLucee, interpSurfMatrixUpper2d[1]);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis2d->getSurfLowerGaussQuadData(1, interpSurfMatrixLower2dLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights2d[1]);
    copyLuceeToEigen(interpSurfMatrixLower2dLucee, interpSurfMatrixLower2d[1]);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates2d[1]);

    // Scale gaussSurfWeights2d[1] for (v,mu) integration
    for (int i = 0; i < gaussSurfWeights2d[1].size(); i++)
      gaussSurfWeights2d[1][i] = gaussSurfWeights2d[1][i]*grid.getDx(3)/grid.getDx(0);

    // Compute and store inverse of mass matrix
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();
    // Compute derivative of basis function evaluated at volume quadrature points
    basisDerivAtVolQuad.resize(2);
    basisDerivAtVolQuad[0] = interpVolMatrix2d*massMatrixInv*gradStiffMatrix[0].transpose();
    basisDerivAtVolQuad[1] = interpVolMatrix2d*massMatrixInv*gradStiffMatrix[1].transpose();
    // Take transpose so same basis function evaluated at diff quad points in the same row
    // Need to use transposeInPlace because otherwise there will be a bug!
    basisDerivAtVolQuad[0].transposeInPlace();
    basisDerivAtVolQuad[1].transposeInPlace();

    // Pre-multiply by inverse mass matrix
    basisDerivAtVolQuad[0] = massMatrixInv*basisDerivAtVolQuad[0];
    basisDerivAtVolQuad[1] = massMatrixInv*basisDerivAtVolQuad[1];

    // Precompute two matrices for surface integrals
    surfIntegralMatrixLower2d.resize(2);
    surfIntegralMatrixUpper2d.resize(2);
    // Pre-multiply by inverse mass matrix
    surfIntegralMatrixLower2d[0] = massMatrixInv*interpSurfMatrixLower2d[0].transpose();
    surfIntegralMatrixUpper2d[0] = massMatrixInv*interpSurfMatrixUpper2d[0].transpose();
    surfIntegralMatrixLower2d[1] = massMatrixInv*interpSurfMatrixLower2d[1].transpose();
    surfIntegralMatrixUpper2d[1] = massMatrixInv*interpSurfMatrixUpper2d[1].transpose();
  }

  Lucee::UpdaterStatus 
  LenardBernsteinDragPara3D2VUpdater::update(double t)
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
    int nlocal2d = nodalBasis2d->getNumNodes();
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
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
    
    // Get alpha. Need to scale by n/T^(3/2) in this class
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    double alpha = resultVector[0];

    // Should be size 4 for linear elements
    Eigen::VectorXd fReduced(nodalStencil.size());
    Eigen::VectorXd fReducedLower(nodalStencil.size());
    Eigen::VectorXd paraQuad(nVolQuad2d);
    Lucee::RowMajorSequencer<5> seq(localRgn);   

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);
      uIn.setPtr(uInPtr, idx[0], idx[1], idx[2]);
      // Set pointers to current location on grid
      fIn.setPtr(fInPtr, idx);
      fOut.setPtr(fOutPtr, idx);
      // Get the coordinates of cell center
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      // At this location, loop over each configuration space grid node
      for (int configNode = 0; configNode < nlocal3d; configNode++)
      {
        // Fill out fReduced at this location
        for (int i = 0; i < fReduced.size(); i++)
          fReduced(i) = fInPtr[configNode + nodalStencil[i]];

        // Compute f at quadrature points
        Eigen::VectorXd fVolQuad = interpVolMatrix2d*fReduced;

        // Compute integrands for gaussian quadrature excluding basis function derivatives
        for (int quadIndex = 0; quadIndex < nVolQuad2d; quadIndex++)
        {
          double paraCoord = cellCentroid[3] + 0.5*grid.getDx(3)*gaussVolOrdinates2d(quadIndex,0);
          paraQuad(quadIndex) = gaussVolWeights2d[quadIndex]*fVolQuad(quadIndex)*
            (paraCoord - uInPtr[configNode]/numDensityInPtr[configNode]);

          // Keep track of max CFL number
          // (from drag in v)
          cfla = std::max(cfla, std::abs(alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))*
              (paraCoord-uInPtr[configNode]/numDensityInPtr[configNode])*dt/grid.getDx(3)));
        }

        // Evaluate integral using gaussian quadrature (represented as matrix-vector multiply)
        Eigen::VectorXd volIntegralResult = basisDerivAtVolQuad[0]*paraQuad;

        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
          fOutPtr[configNode + nodalStencil[nodeIndex]] -= volIntegralResult(nodeIndex);
      }
    }

    // Time-step was too large: return a suggestion with correct time-step
    // Only checking cfl condition at volume quadrature points for now
    if (cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

    // Determine loop bounds for surface integrals
    int ivLower = localRgn.getLower(3);
    // Need one edge outside domain interior
    int ivUpper = localRgn.getUpper(3)+1;
    // Zero flux BC's. Make sure we skip surface integrals on
    // surfaces that lie on a zero-flux boundary
    if (ivLower == globalRgn.getLower(3))
      ivLower = globalRgn.getLower(3)+1;
    if (ivUpper == globalRgn.getUpper(3)+1)
      ivUpper = globalRgn.getUpper(3);

    // Lower/Left, Upper/Right mixed naming convention
    int idxLower[5];
    Eigen::VectorXd surfIntegralFluxes(gaussSurfWeights2d[0].size());

    // Surface integral stage
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      idx[0] = ix;
      idxLower[0] = ix;
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        idx[1] = iy;
        idxLower[1] = iy;
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[2] = iz;
          idxLower[2] = iz;

          temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
          numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);
          uIn.setPtr(uInPtr, idx[0], idx[1], idx[2]);

          // At this location, loop over each configuration space grid node
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // V-Parallel Surface Integrals
            for (int iv = ivLower; iv < ivUpper; iv++)
            {
              idx[3] = iv;
              idxLower[3] = iv-1;
              
              for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
              {
                idx[4] = iMu;
                idxLower[4] = iMu;
                
                // Set index for cell right of interface
                fIn.setPtr(fInPtr, idx);
                // Set index for cell left of interface
                fIn.setPtr(fInPtrLeft, idxLower);

                // Fill out solution for this cell and neighboring cell
                for (int i = 0; i < fReduced.size(); i++)
                {
                  fReduced(i) = fInPtr[configNode + nodalStencil[i]];
                  fReducedLower(i) = fInPtrLeft[configNode + nodalStencil[i]];
                }

                // Evaluate solution at quadrature nodes on same surface
                Eigen::VectorXd fLeftSurfEvals = interpSurfMatrixUpper2d[0]*fReducedLower;
                Eigen::VectorXd fRightSurfEvals = interpSurfMatrixLower2d[0]*fReduced;

                // Need to know center of upper cell to figure out the global velocity coordinate
                grid.setIndex(idx);
                grid.getCentroid(cellCentroid);
                double paraCoord = cellCentroid[3] - 0.5*grid.getDx(3);

                for (int quadIndex = 0; quadIndex < gaussSurfWeights2d[0].size(); quadIndex++)
                {
                  // Compute Lax flux at each surface quadrature point
                  double numFlux;
                  if (paraCoord - uInPtr[configNode]/numDensityInPtr[configNode] > 0.0)
                  {
                    numFlux = (paraCoord - uInPtr[configNode]/numDensityInPtr[configNode])*
                      fRightSurfEvals(quadIndex);
                  }
                  else
                  {
                    numFlux = (paraCoord - uInPtr[configNode]/numDensityInPtr[configNode])*
                      fLeftSurfEvals(quadIndex);
                  }

                  // Store result of weight*(v-uIn)*f at this location in a vector
                  surfIntegralFluxes(quadIndex) = gaussSurfWeights2d[0][quadIndex]*numFlux;
                }
                
                // Compute all surface integrals using a matrix multiply.
                // Each row in result is a basis function times flux projection
                // Then inverse of mass matrix is multiplied to find appropriate increments
                Eigen::VectorXd leftSurfIntegralResult = surfIntegralMatrixUpper2d[0]*surfIntegralFluxes;
                Eigen::VectorXd rightSurfIntegralResult = surfIntegralMatrixLower2d[0]*surfIntegralFluxes;

                // Set index for cell right of interface
                fOut.setPtr(fOutPtr, idx);
                // Set index for cell left of interface
                fOut.setPtr(fOutPtrLeft, idxLower);

                for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                {
                  // Plus sign for left cell since outward normal is in +v direction
                  fOutPtrLeft[configNode+nodalStencil[nodeIndex]] += leftSurfIntegralResult(nodeIndex);
                  // Minus sign for right cell since outward normal is in -v direction
                  fOutPtr[configNode+nodalStencil[nodeIndex]]  -= rightSurfIntegralResult(nodeIndex);
                }
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
                *sqrt(temperatureInPtr[configNode]))*fOutPtr[configNode + nodalStencil[stencilIndex]];
          }
        }
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDragPara3D2VUpdater::declareTypes()
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
  LenardBernsteinDragPara3D2VUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  LenardBernsteinDragPara3D2VUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("LenardBernsteinDragPara3D2VUpdater::evaluateFunction: ");
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
        throw Lucee::Except("LenardBernsteinDragPara3D2VUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

  bool
  LenardBernsteinDragPara3D2VUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

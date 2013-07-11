/**
 * @file	LcLenardBernsteinDiffUpdater.cpp
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinAlgebra.h>
#include <LcLenardBernsteinDiffUpdater.h>


namespace Lucee
{
  using namespace Eigen;
  const char *LenardBernsteinDiffUpdater::id = "LenardBernsteinDiffUpdater2D";

  LenardBernsteinDiffUpdater::LenardBernsteinDiffUpdater()
  {
  }

  LenardBernsteinDiffUpdater::~LenardBernsteinDiffUpdater()
  {
  }

  void
  LenardBernsteinDiffUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("LenardBernsteinDiffUpdater::readInput: Must specify element to use using 'basis'");

    // diffusion coefficient
    alpha = tbl.getNumber("diffusionCoeff");
    // CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number
    // should only increments be computed?
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    // directions to update
    diffDir = 1;
  }

  void
  LenardBernsteinDiffUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();
// allocate space for matrices
    iMat.resize(2);
    lowerMat.resize(2);
    upperMat.resize(2);
    for (unsigned d=0; d<2; ++d)
    {
      iMat[d] = Lucee::Matrix<double>(nlocal, nlocal);
      lowerMat[d] = Lucee::Matrix<double>(nlocal, nlocal);
      upperMat[d] = Lucee::Matrix<double>(nlocal, nlocal);
    }

    Lucee::RowMajorSequencer<2> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[2];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

// get various matrices needed
    nodalBasis->getDiffusionMatrices(iMat, lowerMat, upperMat);

// pre-multiply each of the matrices by inverse matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);
// NOTE: mass matrix is fetched repeatedly as the solve() method
// destroys it during the inversion process
    for (unsigned d=0; d<2; ++d)
    {
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, iMat[d]);
    }
    for (unsigned d=0; d<2; ++d)
    {
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, lowerMat[d]);
    }
    for (unsigned d=0; d<2; ++d)
    {
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, upperMat[d]);
    }

    int numSurfQuadNodes = nodalBasis->getNumSurfGaussNodes();
    Lucee::Matrix<double> interpSurfMatrixLowerLucee(numSurfQuadNodes, nlocal);
    Lucee::Matrix<double> gaussSurfOrdinatesLucee(numSurfQuadNodes, 3);
    gaussSurfWeights = std::vector<double>(numSurfQuadNodes);
    gaussSurfOrdinates = Eigen::MatrixXd(numSurfQuadNodes, 3);
    std::vector<int> lowerSurfNodeNums(nodalBasis->getNumSurfLowerNodes(diffDir));
    Eigen::MatrixXd interpSurfMatrixLower(numSurfQuadNodes, nlocal);
    
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis->getSurfLowerGaussQuadData(diffDir, interpSurfMatrixLowerLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights);
    // Get the nodes on the lower surface. Use this with interpSurfMatrixLucee.
    nodalBasis->getSurfLowerNodeNums(diffDir, lowerSurfNodeNums);
    // Matrix to compute v_t(x)^2 at quadrature locations. Only works in 2-D.
    surfNodeInterpMatrix = Eigen::MatrixXd(numSurfQuadNodes, lowerSurfNodeNums.size());

    copyLuceeToEigen(interpSurfMatrixLowerLucee, interpSurfMatrixLower);

    // Take interpSurfMatrixLower and create a lower dimension interpolation matrix
    for (int nodeIndex = 0; nodeIndex < numSurfQuadNodes; nodeIndex++)
    {
      // At each quadrature node, copy basis function evaluations for
      // those basis functions associated with the nodes on the lower surface
      for (int basisIndex = 0; basisIndex < lowerSurfNodeNums.size(); basisIndex++)
      {
        // Need order of elements in lowerSurfNodeNums to match up with the
        // 1-D data in u(x) for this to work. Not sure how to enforce, but maybe
        // it will just happen to be that way.
        surfNodeInterpMatrix(nodeIndex, basisIndex) = interpSurfMatrixLower(nodeIndex, 
          lowerSurfNodeNums[basisIndex]);
      }
    }
  }

  Lucee::UpdaterStatus
  LenardBernsteinDiffUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    const Lucee::Field<2, double>& inpFld = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<1, double>& inpVtSq = this->getInp<Lucee::Field<1, double> >(1);
    Lucee::Field<2, double>& diffOut = this->getOut<Lucee::Field<2, double> >(0);

    double dt = t-this->getCurrTime();
    if (onlyIncrement)
      diffOut = 0.0;
    else
      diffOut.copy(inpFld);
    
    Lucee::ConstFieldPtr<double> inpFldPtr = inpFld.createConstPtr();
    Lucee::ConstFieldPtr<double> vtSqPtr = inpVtSq.createConstPtr();
    Lucee::FieldPtr<double> diffOutPtr = diffOut.createPtr();

    // check time-step
    double cflm = 1.1*cfl;
    double dxMax2 = grid.getDx(diffDir)*grid.getDx(diffDir);
    double cfla = 0.0;
    double fact = 0.0;

    int idx[2], idxL[2];
// local region to index
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      inpVtSq.setPtr(vtSqPtr, ix);
      Eigen::VectorXd vtSqVals(surfNodeInterpMatrix.cols());
      grid.setIndex(ix, localRgn.getLower(1));
      
      // Copy values of u into uVals vector
      for (int componentIndex = 0; componentIndex < vtSqVals.rows(); componentIndex++)
        vtSqVals(componentIndex) = vtSqPtr[componentIndex];
      
      // Interpolate u to quadrature points on the surface
      // Quadrature coordinates located in gaussSurfOrdinates
      // Each row of uSurfQuad is value of u at the quad location
      Eigen::VectorXd vtSqSurfQuad = surfNodeInterpMatrix*vtSqVals;

      // Integrate vtSqSurf to find its average value
      double vtSqAvg = 0.0;
      for (int quadPoint = 0; quadPoint < vtSqSurfQuad.rows(); quadPoint++)
        vtSqAvg += gaussSurfWeights[quadPoint]*vtSqSurfQuad(quadPoint)/grid.getDx(0);

      // Keep track of maximum cfla
      cfla = std::max(cfla, alpha*vtSqAvg*dt/dxMax2);
      if (cfla > cflm)
        return Lucee::UpdaterStatus(false, dt*cfl/cfla);

      // Update factor
      if (onlyIncrement)
        fact = alpha*vtSqAvg;
      else
        fact = alpha*vtSqAvg*dt;

      for (int iv = localRgn.getLower(1); iv < localRgn.getUpper(1); iv++)
      {
        diffOut.setPtr(diffOutPtr, ix, iv);

        // add in contribution to cell from current cell
        inpFld.setPtr(inpFldPtr, ix, iv);
        matVec(fact, iMat[diffDir], &inpFldPtr[0], 1.0, &diffOutPtr[0]);

        // add in contribution from cells attached to lower/upper faces
        if (iv > globalRgn.getLower(1))
        {
          inpFld.setPtr(inpFldPtr, ix, iv-1); // cell attached to lower face
          matVec(fact, lowerMat[diffDir], &inpFldPtr[0], 1.0, &diffOutPtr[0]);
        }

        if (iv < globalRgn.getUpper(1)-1)
        {
          inpFld.setPtr(inpFldPtr, ix, iv+1); // cell attached to upper face
          matVec(fact, upperMat[diffDir], &inpFldPtr[0], 1.0, &diffOutPtr[0]);
        }
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDiffUpdater::declareTypes()
  {
    // Takes two inputs (f, v_t^2)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void 
  LenardBernsteinDiffUpdater::matVec(double m, const Lucee::Matrix<double>& mat,
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

  void
  LenardBernsteinDiffUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
}

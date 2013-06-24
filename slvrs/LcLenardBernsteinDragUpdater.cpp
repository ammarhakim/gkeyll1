/**
 * @file	LcLenardBernsteinDragUpdater.cpp
 *
 * @brief	Updater to compute the drag term in the L-B collision operator.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcLenardBernsteinDragUpdater.h>
#include <LcStructuredGridBase.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *LenardBernsteinDragUpdater::id = "LenardBernsteinDragUpdater2D";

  LenardBernsteinDragUpdater::LenardBernsteinDragUpdater()
    : UpdaterIfc()
  {
  }

  void 
  LenardBernsteinDragUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("LenardBernsteinDragUpdater::readInput: Must specify element to use using 'basis'");
    
    // diffusion coefficient
    alpha = tbl.getNumber("diffusionCoeff");
    // cfl number
    cfl = tbl.getNumber("cfl");
    // use slightly large max CFL to avoid thrashing around
    cflm = 1.1*cfl;
    // Assume operator will always be in 1 direction
    dragDir = 1;
    
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
    else onlyIncrement = false;
  }

  void 
  LenardBernsteinDragUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<2> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[2];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numVolQuadNodes = nodalBasis->getNumGaussNodes();
    int numSurfQuadNodes = nodalBasis->getNumSurfGaussNodes();

    // Get mass and grad-stiffness matrices and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);
    Lucee::Matrix<double> gradStiffMatrixLucee(nlocal, nlocal);

    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    Eigen::MatrixXd gradStiffMatrix(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    nodalBasis->getGradStiffnessMatrix(dragDir, gradStiffMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpVolMatrixLucee(numVolQuadNodes, nlocal);
    Lucee::Matrix<double> interpSurfMatrixLowerLucee(numSurfQuadNodes, nlocal);
    Lucee::Matrix<double> interpSurfMatrixUpperLucee(numSurfQuadNodes, nlocal);

    Lucee::Matrix<double> gaussVolOrdinatesLucee(numVolQuadNodes, 3);
    Lucee::Matrix<double> gaussSurfOrdinatesLucee(numSurfQuadNodes, 3);

    gaussSurfWeights = std::vector<double>(numSurfQuadNodes);
    gaussVolWeights = std::vector<double>(numVolQuadNodes);
    
    std::vector<int> lowerSurfNodeNums(nodalBasis->getNumSurfLowerNodes(dragDir));

    // Allocate Eigen matrices
    interpVolMatrix = Eigen::MatrixXd(numVolQuadNodes, nlocal);
    interpSurfMatrixLower = Eigen::MatrixXd(numSurfQuadNodes, nlocal);
    interpSurfMatrixUpper = Eigen::MatrixXd(numSurfQuadNodes, nlocal);

    gaussSurfOrdinates = Eigen::MatrixXd(numSurfQuadNodes, 3);
    gaussVolOrdinates = Eigen::MatrixXd(numVolQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpVolMatrixLucee, gaussVolOrdinatesLucee, gaussVolWeights);
    // Get the interpolation matrix for the upper surface quadrature points.
    // Quadrature location and weights will be overwritten but it doesn't matter for these purposes
    nodalBasis->getSurfUpperGaussQuadData(dragDir, interpSurfMatrixUpperLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights);
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis->getSurfLowerGaussQuadData(dragDir, interpSurfMatrixLowerLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights);
    // Get the nodes on the lower surface. Use this with interpSurfMatrixLucee.
    nodalBasis->getSurfLowerNodeNums(dragDir, lowerSurfNodeNums);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix);
    copyLuceeToEigen(interpVolMatrixLucee, interpVolMatrix);
    copyLuceeToEigen(interpSurfMatrixLowerLucee, interpSurfMatrixLower);
    copyLuceeToEigen(interpSurfMatrixUpperLucee, interpSurfMatrixUpper);
    copyLuceeToEigen(gaussSurfOrdinatesLucee, gaussSurfOrdinates);
    copyLuceeToEigen(gaussVolOrdinatesLucee, gaussVolOrdinates);

    // Compute and store inverse of mass matrix
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();
    // Compute derivative of basis function evaluated at volume quadrature points
    basisDerivAtVolQuad = interpVolMatrix*massMatrixInv*gradStiffMatrix.transpose();
    // Take transpose so same basis function evaluated at diff quad points in the same row
    // Need to use transposeInPlace because otherwise there will be a bug!
    basisDerivAtVolQuad.transposeInPlace();
    // Pre-multiply by inverse mass matrix
    basisDerivAtVolQuad = massMatrixInv*basisDerivAtVolQuad;

    surfNodeInterpMatrix = Eigen::MatrixXd(numSurfQuadNodes, lowerSurfNodeNums.size());

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

    // Precompute two matrices for surface integrals
    surfIntegralMatrixLower = massMatrixInv*interpSurfMatrixLower.transpose();
    surfIntegralMatrixUpper = massMatrixInv*interpSurfMatrixUpper.transpose();
  }

  Lucee::UpdaterStatus 
  LenardBernsteinDragUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Distribution function
    const Lucee::Field<2, double>& q = this->getInp<Lucee::Field<2, double> >(0);
    // Drift velocity u(x).
    const Lucee::Field<1, double>& u = this->getInp<Lucee::Field<1, double> >(1);
    // Output distribution function
    Lucee::Field<2, double>& qNew = this->getOut<Lucee::Field<2, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    double dt = t-this->getCurrTime();

    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    double cfla = 0.0; // maximum CFL number used

    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::ConstFieldPtr<double> uPtr = u.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrLeft = qNew.createPtr();

    qNew = 0.0; // use qNew to store increment initially

    int idx[2];
    double cellCentroid[3];

    // Volume integral contribution
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      u.setPtr(uPtr, ix);
      Eigen::VectorXd uVals(surfNodeInterpMatrix.cols());
      
      // Copy values of u into uVals vector
      for (int componentIndex = 0; componentIndex < uVals.rows(); componentIndex++)
        uVals(componentIndex) = uPtr[componentIndex];
      
      // Interpolate u to quadrature points on the surface
      // Quadrature coordinates located in gaussSurfOrdinates
      // Each row of uSurfQuad is value of u at the quad location
      Eigen::VectorXd uSurfQuad = surfNodeInterpMatrix*uVals;

      for (int iv = localRgn.getLower(1); iv < localRgn.getUpper(1); iv++)
      {
        q.setPtr(qPtr, ix, iv);
        qNew.setPtr(qNewPtr, ix, iv);
        grid.setIndex(ix, iv);
        grid.getCentroid(cellCentroid);
        
        // Figure out what f is at each quadrature point
        Eigen::VectorXd fVals(nlocal);
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          fVals(componentIndex) = qPtr[componentIndex];

        // Each row of fVolQuad is f evaluated at a volume quadrature point
        Eigen::VectorXd fVolQuad = interpVolMatrix*fVals;

        // Multiply each element of fVolQuad with (v-u) evaluated at same location
        for (int volNodeIndex = 0; volNodeIndex < fVolQuad.rows(); volNodeIndex++)
        {
          // Loop through gaussSurfOrdinates to find u(x) for this point
          int uMatchedIndex = -1;
          for (int surfNodeIndex = 0; surfNodeIndex < gaussSurfOrdinates.rows(); surfNodeIndex++)
          {
            double numDiff = std::abs(gaussSurfOrdinates(surfNodeIndex,0)-gaussVolOrdinates(volNodeIndex,0));
            if (numDiff < 1.0e-9)
            {
              uMatchedIndex = surfNodeIndex;
              break;
            }
          }

          // Make sure a match for u has been found
          if (uMatchedIndex == -1)
            throw Lucee::Except("LenardBernsteinDragUpdater::update: Unable to find u(x) match in quadrature.");

          double physicalV = cellCentroid[dragDir] + gaussVolOrdinates(volNodeIndex,dragDir)*grid.getDx(dragDir)/2.0;
          fVolQuad(volNodeIndex) = gaussVolWeights[volNodeIndex]*fVolQuad(volNodeIndex)*(physicalV 
            - uSurfQuad(uMatchedIndex));
        }

        // Evaluate integral using gaussian quadrature (represented as matrix-vector multiply)
        Eigen::VectorXd volIntegralResult = basisDerivAtVolQuad*fVolQuad;
        
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          qNewPtr[componentIndex] -= volIntegralResult(componentIndex);
      }
    }

    // Contributions from surface integrals
    // Create sequencer to loop over *each* 1D slice in 'dragDir' direction
    Lucee::RowMajorSequencer<2> seq(localRgn.deflate(dragDir));

    // Lower and upper bounds of 1D slice. (We need to make sure that flux
    // is computed for one edge outside domain interior)
    int sliceLower = localRgn.getLower(dragDir);
    int sliceUpper = localRgn.getUpper(dragDir)+1;

    int idxr[2], idxl[2];
    // Loop over each 1D slice
    // Each seq steps along the 0 direction
    while (seq.step())
    {
      seq.fillWithIndex(idxr);
      seq.fillWithIndex(idxl);

      // Compute u at quadrature points (again)
      u.setPtr(uPtr, idxr[0]);
      Eigen::VectorXd uVals(surfNodeInterpMatrix.cols());

      for (int componentIndex = 0; componentIndex < uVals.rows(); componentIndex++)
        uVals(componentIndex) = uPtr[componentIndex];
      
      // Interpolate u to quadrature points on the surface
      // Quadrature coordinates located in gaussSurfOrdinates
      Eigen::VectorXd uSurfQuad = surfNodeInterpMatrix*uVals;

      // Loop over each edge in slice (in v-direction)
      // Note: I keep the left/right convention, but it makes
      // more sense mentally to go from bottom to top
      for (int i = sliceLower; i < sliceUpper; i++)
      {
        idxr[dragDir] = i; // cell right of edge
        idxl[dragDir] = i-1; // cell left of edge
        
        // Need to know center of right cell to figure out
        // the global velocity coordinate
        grid.setIndex(idxr);
        grid.getCentroid(cellCentroid);
        
        q.setPtr(qPtr, idxr);
        q.setPtr(qPtrl, idxl);
        
        grid.setIndex(idxl);
        double dxL = grid.getDx(dragDir);
        grid.setIndex(idxr);
        double dxR = grid.getDx(dragDir);

        // Copy q data to Eigen vectors for matrix multiplications
        Eigen::VectorXd fLeft(nlocal);
        Eigen::VectorXd fRight(nlocal);
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        {
          fLeft(componentIndex) = qPtrl[componentIndex];
          fRight(componentIndex) = qPtr[componentIndex];
        }

        // Evaluate fLeft and fRight at quadrature nodes on surface
        Eigen::VectorXd fLeftSurfEvals = interpSurfMatrixUpper*fLeft;
        Eigen::VectorXd fRightSurfEvals = interpSurfMatrixLower*fRight;

        // Loop over quadrature points on the edge
        Eigen::VectorXd surfIntegralFluxes(gaussSurfOrdinates.rows());

        // physicalV shouldn't be changing in the loop since all gaussSurfOrdinates should
        // be at a fixed v
        double physicalV = cellCentroid[dragDir] + gaussSurfOrdinates(0,dragDir)*grid.getDx(dragDir)/2.0;
        
        for (int nodeIndex = 0; nodeIndex < gaussSurfOrdinates.rows(); nodeIndex++)
        {
          // Compute Lax flux
          double numFlux = 0.5*(physicalV - uSurfQuad(nodeIndex))*(fLeftSurfEvals(nodeIndex) + fRightSurfEvals(nodeIndex)) -
            0.5*std::abs(physicalV - uSurfQuad(nodeIndex))*(fRightSurfEvals(nodeIndex) - fLeftSurfEvals(nodeIndex));
          // Store result of weight*(v-u)*f at this location in a vector
          surfIntegralFluxes(nodeIndex) = gaussSurfWeights[nodeIndex]*numFlux;

          // Keep track of max CFL number
          cfla = std::max(cfla, alpha*std::abs(physicalV-uSurfQuad(nodeIndex))*dt/grid.getDx(dragDir));
          // Time-step was too large: return a suggestion with correct time-step
          if (cfla > cflm)
            return Lucee::UpdaterStatus(false, dt*cfl/cfla);
        }

        // Compute all surface integrals using a matrix multiply.
        // Each row in result is a basis function times flux projection
        // Then inverse of mass matrix is multiplied to find appropriate increments
        Eigen::VectorXd leftSurfIntegralResult = surfIntegralMatrixUpper*surfIntegralFluxes;
        Eigen::VectorXd rightSurfIntegralResult = surfIntegralMatrixLower*surfIntegralFluxes;

        // Update left cell connected to edge with flux on face
        qNew.setPtr(qNewPtrLeft, idxl);
        // Update right cell connected to edge with flux on face
        qNew.setPtr(qNewPtr, idxr);

        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        {
          // Plus sign for left cell since outward normal is in +v direction
          qNewPtrLeft[componentIndex] += leftSurfIntegralResult(componentIndex);
          // Minus sign for right cell since outward normal is in -v direction
          qNewPtr[componentIndex]  -= rightSurfIntegralResult(componentIndex);
        }
      }
    }
    
    seq = Lucee::RowMajorSequencer<2>(localRgn);
    // Final sweep, update solution with forward Euler step
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      qNew.setPtr(qNewPtr, idx);
      q.setPtr(qPtr, idx);

      // NOTE: If only calculation of increments are requested, the final
      // Euler update is not performed. This means that the multiplication
      // of the DG RHS with dt is not done, something to keep in mind if
      // using the increment in time-dependent update.

      if (onlyIncrement == false)
      {
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          qNewPtr[componentIndex] = qPtr[componentIndex] + dt*alpha*qNewPtr[componentIndex];
      }
      else
      {
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          qNewPtr[componentIndex] = alpha*qNewPtr[componentIndex];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDragUpdater::declareTypes()
  {
    // takes two inputs (fOld, u) 
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output (fNew)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  LenardBernsteinDragUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

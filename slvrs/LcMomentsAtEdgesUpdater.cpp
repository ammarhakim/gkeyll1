/**
 * @file	LcMomentsAtEdgesUpdater.cpp
 *
 * @brief	Computes first and third moments of the distribution function at both edges.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMomentsAtEdgesUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  const char *MomentsAtEdgesUpdater::id = "MomentsAtEdgesUpdater";

  MomentsAtEdgesUpdater::MomentsAtEdgesUpdater()
  {
  }

  void
  MomentsAtEdgesUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("MomentsAtEdgesUpdater::readInput: Must specify element to use using 'basis'");
  }

  void
  MomentsAtEdgesUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned numNodes = nodalBasis->getNumNodes();
    std::vector<unsigned> yRef(numNodes), xRef(numNodes);

    // reflection mapping
    rotMap.resize(numNodes);
    nodalBasis->getUpperReflectingBcMapping(0, yRef);
    nodalBasis->getLowerReflectingBcMapping(1, xRef);
    for (unsigned i=0; i<numNodes; ++i)
      rotMap[i] = xRef[yRef[i]];

    // Figure out what nodes are on the right edge
    int numEdgeQuadNodes = nodalBasis->getNumSurfGaussNodes();
    Lucee::Matrix<double> interpEdgeMatrixLucee(numEdgeQuadNodes, numNodes);
    Lucee::Matrix<double> gaussEdgeOrdinatesLucee(numEdgeQuadNodes, 3);
    // Allocate Eigen matrices
    Eigen::MatrixXd interpEdgeMatrix(numEdgeQuadNodes, numNodes);
    gaussEdgeOrdinates = Eigen::MatrixXd(numEdgeQuadNodes, 3);
    gaussEdgeWeights = std::vector<double>(numEdgeQuadNodes);
    // Get the interpolation matrix for the right edge quadrature points.
    nodalBasis->getSurfUpperGaussQuadData(0, interpEdgeMatrixLucee, gaussEdgeOrdinatesLucee,
      gaussEdgeWeights);

    rightEdgeNodeNums = std::vector<int>(nodalBasis->getNumSurfUpperNodes(0));
    leftEdgeNodeNums = std::vector<int>(nodalBasis->getNumSurfLowerNodes(0));
    nodalBasis->getSurfUpperNodeNums(0, rightEdgeNodeNums);
    nodalBasis->getSurfLowerNodeNums(0, leftEdgeNodeNums);

    // Copy matrices to eigen objects
    copyLuceeToEigen(interpEdgeMatrixLucee, interpEdgeMatrix);
    copyLuceeToEigen(gaussEdgeOrdinatesLucee, gaussEdgeOrdinates);

    edgeNodeInterpMatrix = Eigen::MatrixXd(numEdgeQuadNodes, rightEdgeNodeNums.size());

    // Create a lower dimension (1-D) interpolation matrix
    for (int nodeIndex = 0; nodeIndex < numEdgeQuadNodes; nodeIndex++)
    {
      // At each quadrature node, copy basis function evaluations for
      // those basis functions associated with the nodes on the (right) edge
      for (int basisIndex = 0; basisIndex < rightEdgeNodeNums.size(); basisIndex++)
      {
        edgeNodeInterpMatrix(nodeIndex, basisIndex) = interpEdgeMatrix(nodeIndex, 
          rightEdgeNodeNums[basisIndex]);
      }
    }

    // Derivative stuff
    Lucee::Matrix<double> massMatrixLucee(numNodes, numNodes);
    Eigen::MatrixXd massMatrix(numNodes, numNodes);
    nodalBasis->getMassMatrix(massMatrixLucee);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    
    Lucee::Matrix<double> gradStiffMatrixLucee(numNodes, numNodes);
    Eigen::MatrixXd gradStiffMatrix(numNodes, numNodes);
    nodalBasis->getGradStiffnessMatrix(1, gradStiffMatrixLucee);
    copyLuceeToEigen(gradStiffMatrixLucee, gradStiffMatrix);

    derivativeMatrix = massMatrix.inverse()*gradStiffMatrix.transpose();
  }

  Lucee::UpdaterStatus
  MomentsAtEdgesUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Distribution function
    const Lucee::Field<2, double>& distfIn = this->getInp<Lucee::Field<2, double> >(0);
    // Hamiltonian
    const Lucee::Field<2, double>& hamilIn = this->getInp<Lucee::Field<2, double> >(1);
    // Output velocity moments (1 and 3) vs time
    Lucee::DynVector<double>& velocityMoments = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();
    
    Lucee::ConstFieldPtr<double> sknPtr = distfIn.createConstPtr(); // for skin-cell
    Lucee::ConstFieldPtr<double> gstPtr = distfIn.createConstPtr(); // for ghost-cell
    Lucee::ConstFieldPtr<double> sknHamilPtr = hamilIn.createConstPtr(); // for skin-cell
    Lucee::ConstFieldPtr<double> gstHamilPtr = hamilIn.createConstPtr(); // for ghost-cell

    unsigned numNodes = nodalBasis->getNumNodes();

    double mom0AlongRightEdge = 0.0;
    double mom0AlongLeftEdge = 0.0;
    double mom1AlongRightEdge = 0.0;
    double mom1AlongLeftEdge = 0.0;
    double mom2AlongRightEdge = 0.0;
    double mom2AlongLeftEdge = 0.0;
    double mom3AlongRightEdge = 0.0;
    double mom3AlongLeftEdge = 0.0;

    double cellCentroid[3];
    int idx[2];

    int ix = globalRgn.getUpper(0)-1; // right skin cell x index
    
    // Integrate along skin cell (right edge)
    for (int js = globalRgn.getUpper(1)-1; js >= 0; js--)
    {
      idx[0] = ix;
      idx[1] = js;
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      // Don't want skin cell contributions for v < 0.0 (upwinding)
      if (cellCentroid[1] < 0.0)
        break;

      distfIn.setPtr(sknPtr, ix, js);
      hamilIn.setPtr(sknHamilPtr, ix, js);

      Eigen::VectorXd hamilAtNodes(numNodes);
      for (int i = 0; i < numNodes; i++)
        hamilAtNodes(i) = sknHamilPtr[i];
      Eigen::VectorXd hamilDerivAtNodes = derivativeMatrix*hamilAtNodes;

      // Copy nodes on right edge into a vector
      Eigen::VectorXd edgeHamilDerivData(rightEdgeNodeNums.size());
      Eigen::VectorXd edgeHamilData(rightEdgeNodeNums.size());
      Eigen::VectorXd edgeData(rightEdgeNodeNums.size());

      // Copy nodal values of hamiltonian, hamiltonian derivative, and dist f into a vector
      for (int edgeNodeIndex = 0; edgeNodeIndex < edgeHamilData.rows(); edgeNodeIndex++)
      {
        edgeHamilData(edgeNodeIndex)      = hamilAtNodes(rightEdgeNodeNums[edgeNodeIndex]);
        edgeHamilDerivData(edgeNodeIndex) = hamilDerivAtNodes(rightEdgeNodeNums[edgeNodeIndex]);
        edgeData(edgeNodeIndex) = sknPtr[rightEdgeNodeNums[edgeNodeIndex]];
      }

      // Interpolate nodal data to quadrature points on the edge
      Eigen::VectorXd edgeHamilQuadData = edgeNodeInterpMatrix*edgeHamilData;
      Eigen::VectorXd edgeHamilDerivQuadData = edgeNodeInterpMatrix*edgeHamilDerivData;
      Eigen::VectorXd edgeQuadData = edgeNodeInterpMatrix*edgeData;

      // Integrate v*f over entire cell using gaussian quadrature
      double mom0InEntireCell = 0.0;
      double mom1InEntireCell = 0.0;
      double mom2InEntireCell = 0.0;
      double mom3InEntireCell = 0.0;

      for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
      {
        double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
        
        mom0InEntireCell += gaussEdgeWeights[quadNodeIndex]*edgeQuadData(quadNodeIndex);

        mom1InEntireCell += gaussEdgeWeights[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

        mom2InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeQuadData[quadNodeIndex];

        mom3InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeHamilDerivQuadData[quadNodeIndex]*edgeQuadData[quadNodeIndex];
      }

      mom0AlongRightEdge += mom0InEntireCell;
      mom1AlongRightEdge += mom1InEntireCell;
      mom2AlongRightEdge += mom2InEntireCell;
      mom3AlongRightEdge += mom3InEntireCell;
    }

    // Integrate along ghost cell (right edge)
    for (int jg = globalRgn.getLower(1); jg < globalRgn.getUpper(1); jg++)
    {
      idx[0] = ix+1;
      idx[1] = jg;
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      // Don't want ghost cell contributions for v > 0.0 (upwinding)
      if (cellCentroid[1] > 0.0)
        break;

      distfIn.setPtr(gstPtr, ix+1, jg);
      hamilIn.setPtr(gstHamilPtr, ix, jg);

      Eigen::VectorXd hamilAtNodes(numNodes);
      for (int i = 0; i < numNodes; i++)
        hamilAtNodes(i) = gstHamilPtr[i];
      Eigen::VectorXd hamilDerivAtNodes = derivativeMatrix*hamilAtNodes;

      // Copy nodes on left edge into a vector
      Eigen::VectorXd edgeHamilDerivData(leftEdgeNodeNums.size());
      Eigen::VectorXd edgeHamilData(leftEdgeNodeNums.size());
      Eigen::VectorXd edgeData(leftEdgeNodeNums.size());

      // Copy nodal values of hamiltonian, hamiltonian derivative, and dist f into a vector
      for (int edgeNodeIndex = 0; edgeNodeIndex < edgeHamilData.rows(); edgeNodeIndex++)
      {
        edgeHamilData(edgeNodeIndex)      = hamilAtNodes(rightEdgeNodeNums[edgeNodeIndex]);
        edgeHamilDerivData(edgeNodeIndex) = hamilDerivAtNodes(rightEdgeNodeNums[edgeNodeIndex]);
        edgeData(edgeNodeIndex) = gstPtr[leftEdgeNodeNums[edgeNodeIndex]];
      }

      // Interpolate nodal data to quadrature points on the edge
      Eigen::VectorXd edgeHamilQuadData = edgeNodeInterpMatrix*edgeHamilData;
      Eigen::VectorXd edgeHamilDerivQuadData = edgeNodeInterpMatrix*edgeHamilDerivData;
      Eigen::VectorXd edgeQuadData = edgeNodeInterpMatrix*edgeData;

      // Integrate v*f over entire cell using gaussian quadrature
      double mom0InEntireCell = 0.0;
      double mom1InEntireCell = 0.0;
      double mom2InEntireCell = 0.0;
      double mom3InEntireCell = 0.0;

      for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
      {
        double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
        
        mom0InEntireCell += gaussEdgeWeights[quadNodeIndex]*edgeQuadData(quadNodeIndex);
        
        mom1InEntireCell += gaussEdgeWeights[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

        mom2InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeQuadData[quadNodeIndex];

        mom3InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeHamilDerivQuadData[quadNodeIndex]*edgeQuadData[quadNodeIndex];
      }

      mom0AlongRightEdge += mom0InEntireCell;
      mom1AlongRightEdge += mom1InEntireCell;
      mom2AlongRightEdge += mom2InEntireCell;
      mom3AlongRightEdge += mom3InEntireCell;
    }

    ix = globalRgn.getLower(0); // left skin cell x index

    // Integrate along skin cell (left edge)
    for (int js = 0; js < globalRgn.getUpper(1); js++)
    {
      idx[0] = ix;
      idx[1] = js;
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      // Don't want skin cell contributions for v > 0.0 (upwinding)
      if (cellCentroid[1] > 0.0)
        break;

      distfIn.setPtr(sknPtr, ix, js);
      hamilIn.setPtr(sknHamilPtr, ix, js);

      Eigen::VectorXd hamilAtNodes(numNodes);
      for (int i = 0; i < numNodes; i++)
        hamilAtNodes(i) = sknHamilPtr[i];
      Eigen::VectorXd hamilDerivAtNodes = derivativeMatrix*hamilAtNodes;

      // Copy nodes on left edge into a vector
      Eigen::VectorXd edgeHamilDerivData(leftEdgeNodeNums.size());
      Eigen::VectorXd edgeHamilData(leftEdgeNodeNums.size());
      Eigen::VectorXd edgeData(leftEdgeNodeNums.size());

      // Copy nodal values of hamiltonian, hamiltonian derivative, and dist f into a vector
      for (int edgeNodeIndex = 0; edgeNodeIndex < edgeHamilData.rows(); edgeNodeIndex++)
      {
        edgeHamilData(edgeNodeIndex)      = hamilAtNodes(leftEdgeNodeNums[edgeNodeIndex]);
        edgeHamilDerivData(edgeNodeIndex) = hamilDerivAtNodes(leftEdgeNodeNums[edgeNodeIndex]);
        edgeData(edgeNodeIndex) = sknPtr[leftEdgeNodeNums[edgeNodeIndex]];
      }

      // Interpolate nodal data to quadrature points on the edge
      Eigen::VectorXd edgeHamilQuadData = edgeNodeInterpMatrix*edgeHamilData;
      Eigen::VectorXd edgeHamilDerivQuadData = edgeNodeInterpMatrix*edgeHamilDerivData;
      Eigen::VectorXd edgeQuadData = edgeNodeInterpMatrix*edgeData;

      // Integrate v*f over entire cell using gaussian quadrature
      double mom0InEntireCell = 0.0;
      double mom1InEntireCell = 0.0;
      double mom2InEntireCell = 0.0;
      double mom3InEntireCell = 0.0;

      for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
      {
        double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
        
        mom0InEntireCell += gaussEdgeWeights[quadNodeIndex]*edgeQuadData(quadNodeIndex);
        
        mom1InEntireCell += gaussEdgeWeights[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

        mom2InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeQuadData[quadNodeIndex];

        mom3InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeHamilDerivQuadData[quadNodeIndex]*edgeQuadData[quadNodeIndex];
      }

      mom0AlongLeftEdge += mom0InEntireCell;
      mom1AlongLeftEdge += mom1InEntireCell;
      mom2AlongLeftEdge += mom2InEntireCell;
      mom3AlongLeftEdge += mom3InEntireCell;
    }

    // Integrate along ghost cell (left edge)
    for (int jg = globalRgn.getUpper(1)-1; jg > 0; jg--)
    {
      idx[0] = ix-1;
      idx[1] = jg;
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      // Don't want ghost cell contributions for v < 0.0 (upwinding)
      if (cellCentroid[1] < 0.0)
        break;

      distfIn.setPtr(gstPtr, ix-1, jg);
      hamilIn.setPtr(gstHamilPtr, ix, jg);

      Eigen::VectorXd hamilAtNodes(numNodes);
      for (int i = 0; i < numNodes; i++)
        hamilAtNodes(i) = gstHamilPtr[i];
      Eigen::VectorXd hamilDerivAtNodes = derivativeMatrix*hamilAtNodes;

      // Copy nodes on right edge into a vector
      Eigen::VectorXd edgeHamilDerivData(rightEdgeNodeNums.size());
      Eigen::VectorXd edgeHamilData(rightEdgeNodeNums.size());
      Eigen::VectorXd edgeData(rightEdgeNodeNums.size());

      // Copy nodal values of hamiltonian, hamiltonian derivative, and dist f into a vector
      for (int edgeNodeIndex = 0; edgeNodeIndex < edgeHamilData.rows(); edgeNodeIndex++)
      {
        edgeHamilData(edgeNodeIndex)      = hamilAtNodes(leftEdgeNodeNums[edgeNodeIndex]);
        edgeHamilDerivData(edgeNodeIndex) = hamilDerivAtNodes(leftEdgeNodeNums[edgeNodeIndex]);
        edgeData(edgeNodeIndex) = gstPtr[rightEdgeNodeNums[edgeNodeIndex]];
      }

      // Interpolate nodal data to quadrature points on the edge
      Eigen::VectorXd edgeHamilQuadData = edgeNodeInterpMatrix*edgeHamilData;
      Eigen::VectorXd edgeHamilDerivQuadData = edgeNodeInterpMatrix*edgeHamilDerivData;
      Eigen::VectorXd edgeQuadData = edgeNodeInterpMatrix*edgeData;

      // Integrate v*f over entire cell using gaussian quadrature
      double mom0InEntireCell = 0.0;
      double mom1InEntireCell = 0.0;
      double mom2InEntireCell = 0.0;
      double mom3InEntireCell = 0.0;

      for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
      {
        double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
        
        mom0InEntireCell += gaussEdgeWeights[quadNodeIndex]*edgeQuadData(quadNodeIndex);
        
        mom1InEntireCell += gaussEdgeWeights[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

        mom2InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeQuadData[quadNodeIndex];

        mom3InEntireCell += gaussEdgeWeights[quadNodeIndex]*2*edgeHamilQuadData[quadNodeIndex]*
          edgeHamilDerivQuadData[quadNodeIndex]*edgeQuadData[quadNodeIndex];
      }

      mom0AlongLeftEdge += mom0InEntireCell;
      mom1AlongLeftEdge += mom1InEntireCell;
      mom2AlongLeftEdge += mom2InEntireCell;
      mom3AlongLeftEdge += mom3InEntireCell;
    }

    std::vector<double> data(8);
    
    data[0] = mom0AlongLeftEdge;
    data[1] = mom1AlongLeftEdge;
    data[2] = mom2AlongLeftEdge;
    data[3] = mom3AlongLeftEdge;
    data[4] = mom0AlongRightEdge;
    data[5] = mom1AlongRightEdge;
    data[6] = mom2AlongRightEdge;
    data[7] = mom3AlongRightEdge;

    // Put data into the DynVector
    velocityMoments.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  MomentsAtEdgesUpdater::declareTypes()
  {
    // A distribution function
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    // A hamiltonian
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    // Moments 0-3 at left and right edges
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  MomentsAtEdgesUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

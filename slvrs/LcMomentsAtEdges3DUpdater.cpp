/**
 * @file	LcMomentsAtEdges3DUpdater.cpp
 *
 * @brief	Computes several parallel velocity moments of the distribution function at both edges.
 * Used for 3D SOL problem (Adiabatic Electrons).
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMomentsAtEdges3DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  const char *MomentsAtEdges3DUpdater::id = "MomentsAtEdges3DUpdater";

  MomentsAtEdges3DUpdater::MomentsAtEdges3DUpdater()
  {
  }

  void
  MomentsAtEdges3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("MomentsAtEdges3DUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  MomentsAtEdges3DUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned numNodes = nodalBasis->getNumNodes();

    // Figure out what nodes are on the right edge
    int numEdgeQuadNodes = nodalBasis->getNumSurfGaussNodes();
    Lucee::Matrix<double> interpEdgeMatrixLucee(numEdgeQuadNodes, numNodes);
    Lucee::Matrix<double> gaussEdgeOrdinatesLucee(numEdgeQuadNodes, 3);
    // Allocate Eigen matrices
    interpEdgeMatrixLower = Eigen::MatrixXd(numEdgeQuadNodes, numNodes);
    interpEdgeMatrixUpper = Eigen::MatrixXd(numEdgeQuadNodes, numNodes);
    
    gaussEdgeOrdinatesLower = Eigen::MatrixXd(numEdgeQuadNodes, 3);
    gaussEdgeWeightsLower = std::vector<double>(numEdgeQuadNodes);
    gaussEdgeOrdinatesUpper = Eigen::MatrixXd(numEdgeQuadNodes, 3);
    gaussEdgeWeightsUpper = std::vector<double>(numEdgeQuadNodes);

    // Get the interpolation matrix for the right edge of the domain
    nodalBasis->getSurfUpperGaussQuadData(0, interpEdgeMatrixLucee, gaussEdgeOrdinatesLucee,
      gaussEdgeWeightsUpper);
    // Copy matrices to eigen objects
    copyLuceeToEigen(interpEdgeMatrixLucee, interpEdgeMatrixUpper);
    copyLuceeToEigen(gaussEdgeOrdinatesLucee, gaussEdgeOrdinatesUpper);

    // Do the same for the left edge of the domain
    nodalBasis->getSurfLowerGaussQuadData(0, interpEdgeMatrixLucee, gaussEdgeOrdinatesLucee,
      gaussEdgeWeightsLower);
    // Copy matrices to eigen objects
    copyLuceeToEigen(interpEdgeMatrixLucee, interpEdgeMatrixLower);
    copyLuceeToEigen(gaussEdgeOrdinatesLucee, gaussEdgeOrdinatesLower);

    // Derivative stuff (all unused)
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
  MomentsAtEdges3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Distribution function
    const Lucee::Field<3, double>& distfIn = this->getInp<Lucee::Field<3, double> >(0);
    // Hamiltonian
    const Lucee::Field<3, double>& hamilIn = this->getInp<Lucee::Field<3, double> >(1);
    // Output velocity moments (0, 1, 2, 3) vs time
    Lucee::DynVector<double>& velocityMoments = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    
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
    int idx[3];

    int ix = globalRgn.getUpper(0)-1; // right skin cell x index
    
    // Integrate along skin cell (right edge)
    for (int js = globalRgn.getUpper(1)-1; js >= 0; js--)
    {
      for (int iMu = globalRgn.getLower(2); iMu < globalRgn.getUpper(2); iMu++)
      {
        idx[0] = ix;
        idx[1] = js;
        idx[2] = iMu;

        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't want skin cell contributions for vPara < 0.0 (upwinding)
        if (cellCentroid[1] < 0.0)
          break;

        distfIn.setPtr(sknPtr, ix, js, iMu);
        hamilIn.setPtr(sknHamilPtr, ix, js, iMu);

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = sknPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixUpper*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double mom0InEntireCell = 0.0;
        double mom1InEntireCell = 0.0;
        double mom2InEntireCell = 0.0;
        double mom3InEntireCell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesUpper(quadNodeIndex,1)*grid.getDx(1)/2.0;
          
          mom0InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*edgeQuadData(quadNodeIndex);

          mom1InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

          mom2InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalV*edgeQuadData(quadNodeIndex);

          mom3InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData(quadNodeIndex);
        }
        
        mom0AlongRightEdge += mom0InEntireCell;
        mom1AlongRightEdge += mom1InEntireCell;
        mom2AlongRightEdge += mom2InEntireCell;
        mom3AlongRightEdge += mom3InEntireCell;
      }
    }

    
    // Integrate along ghost cell (right edge)
    for (int jg = globalRgn.getLower(1); jg < globalRgn.getUpper(1); jg++)
    {
      for (int iMu = globalRgn.getLower(2); iMu < globalRgn.getUpper(2); iMu++)
      {
        idx[0] = ix+1;
        idx[1] = jg;
        idx[2] = iMu;

        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't want ghost cell contributions for v > 0.0 (upwinding)
        if (cellCentroid[1] > 0.0)
          break;

        distfIn.setPtr(gstPtr, ix+1, jg, iMu);
        hamilIn.setPtr(gstHamilPtr, ix, jg, iMu);

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = gstPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixLower*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double mom0InEntireCell = 0.0;
        double mom1InEntireCell = 0.0;
        double mom2InEntireCell = 0.0;
        double mom3InEntireCell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesLower(quadNodeIndex,1)*grid.getDx(1)/2.0;
          
          mom0InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*edgeQuadData(quadNodeIndex);
          
          mom1InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

          mom2InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalV*edgeQuadData[quadNodeIndex];

          mom3InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData[quadNodeIndex];
        }

        mom0AlongRightEdge += mom0InEntireCell;
        mom1AlongRightEdge += mom1InEntireCell;
        mom2AlongRightEdge += mom2InEntireCell;
        mom3AlongRightEdge += mom3InEntireCell;
      }
    }

    ix = globalRgn.getLower(0); // left skin cell x index

    // Integrate along skin cell (left edge)
    for (int js = 0; js < globalRgn.getUpper(1); js++)
    {
      for (int iMu = globalRgn.getLower(2); iMu < globalRgn.getUpper(2); iMu++)
      {
        idx[0] = ix;
        idx[1] = js;
        idx[2] = iMu;

        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't want skin cell contributions for v > 0.0 (upwinding)
        if (cellCentroid[1] > 0.0)
          break;

        distfIn.setPtr(sknPtr, ix, js, iMu);
        hamilIn.setPtr(sknHamilPtr, ix, js, iMu);

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = sknPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixLower*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double mom0InEntireCell = 0.0;
        double mom1InEntireCell = 0.0;
        double mom2InEntireCell = 0.0;
        double mom3InEntireCell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesLower(quadNodeIndex,1)*grid.getDx(1)/2.0;
          
          mom0InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*edgeQuadData(quadNodeIndex);
          
          mom1InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

          mom2InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalV*edgeQuadData[quadNodeIndex];

          mom3InEntireCell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData[quadNodeIndex];
        }

        mom0AlongLeftEdge += mom0InEntireCell;
        mom1AlongLeftEdge += mom1InEntireCell;
        mom2AlongLeftEdge += mom2InEntireCell;
        mom3AlongLeftEdge += mom3InEntireCell;
      }
    }

    // Integrate along ghost cell (left edge)
    for (int jg = globalRgn.getUpper(1)-1; jg > 0; jg--)
    {
      for (int iMu = globalRgn.getLower(2); iMu < globalRgn.getUpper(2); iMu++)
      {
        idx[0] = ix-1;
        idx[1] = jg;
        idx[2] = iMu;

        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't want ghost cell contributions for v < 0.0 (upwinding)
        if (cellCentroid[1] < 0.0)
          break;

        distfIn.setPtr(gstPtr, ix-1, jg, iMu);
        hamilIn.setPtr(gstHamilPtr, ix, jg, iMu);

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = gstPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixUpper*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double mom0InEntireCell = 0.0;
        double mom1InEntireCell = 0.0;
        double mom2InEntireCell = 0.0;
        double mom3InEntireCell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesUpper(quadNodeIndex,1)*grid.getDx(1)/2.0;
          
          mom0InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*edgeQuadData(quadNodeIndex);
          
          mom1InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);

          mom2InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalV*edgeQuadData[quadNodeIndex];

          mom3InEntireCell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData[quadNodeIndex];
        }

        mom0AlongLeftEdge += mom0InEntireCell;
        mom1AlongLeftEdge += mom1InEntireCell;
        mom2AlongLeftEdge += mom2InEntireCell;
        mom3AlongLeftEdge += mom3InEntireCell;
      }
    }

    std::vector<double> data(8);

    data[0] = scaleFactor*mom0AlongLeftEdge;
    data[1] = scaleFactor*mom1AlongLeftEdge;
    data[2] = scaleFactor*mom2AlongLeftEdge;
    data[3] = scaleFactor*mom3AlongLeftEdge;
    data[4] = scaleFactor*mom0AlongRightEdge;
    // Used in heat flux calculation:
    data[5] = scaleFactor*mom1AlongRightEdge;
    data[6] = scaleFactor*mom2AlongRightEdge;
    // Used in heat flux calculation:
    data[7] = scaleFactor*mom3AlongRightEdge;

    // Put data into the DynVector
    velocityMoments.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  MomentsAtEdges3DUpdater::declareTypes()
  {
    // A distribution function
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // A hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Moments 0-3 at left and right edges
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  MomentsAtEdges3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

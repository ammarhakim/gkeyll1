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
    // Output velocity moments (0, 1, 2, 3) vs time
    Lucee::DynVector<double>& velocityMoments = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    
    Lucee::ConstFieldPtr<double> sknPtr = distfIn.createConstPtr(); // for skin-cell
    Lucee::ConstFieldPtr<double> gstPtr = distfIn.createConstPtr(); // for ghost-cell

    unsigned numNodes = nodalBasis->getNumNodes();

    //double mom0AlongLeftEdge = 0.0;
    //double mom1AlongLeftEdge = 0.0;
    //double mom2AlongLeftEdge = 0.0;
    //double mom3AlongLeftEdge = 0.0;
    
    //<vPara*f>
    double vParaMom_r = 0.0;
    //<vPara^3*f>
    double vPara3Mom_r = 0.0;
    //<vPara*mu*f>
    double vParaMuMom_r = 0.0;

    double cellCentroid[3];
    int idx[3];

    int ix = globalRgn.getUpper(0)-1; // right skin cell x index
    
    // Integrate along skin cell (right edge)
    for (int js = globalRgn.getLower(1); js < globalRgn.getUpper(1); js++)
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

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = sknPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixUpper*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double vParaMom_cell = 0.0;
        double vPara3Mom_cell = 0.0;
        double vParaMuMom_cell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesUpper(quadNodeIndex,1)*grid.getDx(1)/2.0;
          double physicalMu = cellCentroid[2] + gaussEdgeOrdinatesUpper(quadNodeIndex,2)*grid.getDx(2)/2.0;

          vParaMom_cell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);
          vPara3Mom_cell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData(quadNodeIndex);
          vParaMuMom_cell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalMu*edgeQuadData(quadNodeIndex);
        }
        
        vParaMom_r += vParaMom_cell;
        vPara3Mom_r += vPara3Mom_cell;
        vParaMuMom_r += vParaMuMom_cell;
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

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = gstPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixLower*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double vParaMom_cell = 0.0;
        double vPara3Mom_cell = 0.0;
        double vParaMuMom_cell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesLower(quadNodeIndex,1)*grid.getDx(1)/2.0;
          double physicalMu = cellCentroid[2] + gaussEdgeOrdinatesLower(quadNodeIndex,2)*grid.getDx(2)/2.0;
          
          vParaMom_cell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);
          vPara3Mom_cell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData(quadNodeIndex);
          vParaMuMom_cell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalMu*edgeQuadData(quadNodeIndex);
        }
        
        vParaMom_r += vParaMom_cell;
        vPara3Mom_r += vPara3Mom_cell;
        vParaMuMom_r += vParaMuMom_cell;
      }
    }

    //<vPara*f>
    double vParaMom_l = 0.0;
    //<vPara^3*f>
    double vPara3Mom_l = 0.0;
    //<vPara*mu*f>
    double vParaMuMom_l = 0.0;

    ix = globalRgn.getLower(0); // left skin cell x index
    
    // Integrate along skin cell (left edge)
    for (int js = globalRgn.getLower(1); js < globalRgn.getUpper(1); js++)
    {
      for (int iMu = globalRgn.getLower(2); iMu < globalRgn.getUpper(2); iMu++)
      {
        idx[0] = ix;
        idx[1] = js;
        idx[2] = iMu;

        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't want skin cell contributions for vPara > 0.0 (upwinding)
        if (cellCentroid[1] > 0.0)
          break;

        distfIn.setPtr(sknPtr, ix, js, iMu);

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = sknPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixLower*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double vParaMom_cell = 0.0;
        double vPara3Mom_cell = 0.0;
        double vParaMuMom_cell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesLower(quadNodeIndex,1)*grid.getDx(1)/2.0;
          double physicalMu = cellCentroid[2] + gaussEdgeOrdinatesLower(quadNodeIndex,2)*grid.getDx(2)/2.0;

          vParaMom_cell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);
          vPara3Mom_cell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData(quadNodeIndex);
          vParaMuMom_cell += gaussEdgeWeightsLower[quadNodeIndex]*physicalV*physicalMu*edgeQuadData(quadNodeIndex);
        }
        
        vParaMom_l += vParaMom_cell;
        vPara3Mom_l += vPara3Mom_cell;
        vParaMuMom_l += vParaMuMom_cell;
      }
    }

    // Integrate along ghost cell (left edge)
    for (int jg = globalRgn.getLower(1); jg < globalRgn.getUpper(1); jg++)
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

        distfIn.setPtr(gstPtr, ix+1, jg, iMu);

        // Interpolate f to surface quadrature points
        Eigen::VectorXd fAtNodes(numNodes);
        for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
          fAtNodes(nodeIndex) = gstPtr[nodeIndex];
        Eigen::VectorXd edgeQuadData = interpEdgeMatrixUpper*fAtNodes;

        // Integrate v*f over entire cell using gaussian quadrature
        double vParaMom_cell = 0.0;
        double vPara3Mom_cell = 0.0;
        double vParaMuMom_cell = 0.0;

        for (int quadNodeIndex = 0; quadNodeIndex < edgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinatesUpper(quadNodeIndex,1)*grid.getDx(1)/2.0;
          double physicalMu = cellCentroid[2] + gaussEdgeOrdinatesUpper(quadNodeIndex,2)*grid.getDx(2)/2.0;
          
          vParaMom_cell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*edgeQuadData(quadNodeIndex);
          vPara3Mom_cell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalV*physicalV*edgeQuadData(quadNodeIndex);
          vParaMuMom_cell += gaussEdgeWeightsUpper[quadNodeIndex]*physicalV*physicalMu*edgeQuadData(quadNodeIndex);
        }
        
        vParaMom_l += vParaMom_cell;
        vPara3Mom_l += vPara3Mom_cell;
        vParaMuMom_l += vParaMuMom_cell;
      }
    }



    std::vector<double> data(6);
    data[0] = scaleFactor*vParaMom_r;
    data[1] = scaleFactor*vPara3Mom_r;
    data[2] = scaleFactor*vParaMuMom_r;
    data[3] = scaleFactor*vParaMom_l;
    data[4] = scaleFactor*vPara3Mom_l;
    data[5] = scaleFactor*vParaMuMom_l;

    // Put data into the DynVector
    velocityMoments.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  MomentsAtEdges3DUpdater::declareTypes()
  {
    // A distribution function
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Moments 0-3 at right and left edges
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

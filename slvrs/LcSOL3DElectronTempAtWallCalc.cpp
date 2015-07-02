/**
 * @file	LcSOL3DElectronTempAtWallCalc.cpp
 *
 * @brief	Used to compute electron temperature at wall when sheath boundary
 * conditions are used by taking into account potential drop and cutoff velocity
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOL3DElectronTempAtWallCalc.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  const char *SOL3DElectronTempAtWallCalc::id = "SOL3DElectronTempAtWallCalc";

  SOL3DElectronTempAtWallCalc::SOL3DElectronTempAtWallCalc()
  {
  }

  void
  SOL3DElectronTempAtWallCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOL3DElectronTempAtWallCalc::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("elcMass"))
      elcMass = tbl.getNumber("elcMass");
    else
      throw Lucee::Except("SOL3DElectronTempAtWallCalc::readInput: Must specify electron mass using 'elcMass'");

    if (tbl.hasNumber("eV"))
      eV = tbl.getNumber("eV");
    else
      throw Lucee::Except("SOL3DElectronTempAtWallCalc::readInput: Must specify units of eV using 'eV'");

    if (tbl.hasNumber("B0"))
      B0 = tbl.getNumber("B0");
    else
      throw Lucee::Except("SOL3DElectronTempAtWallCalc::readInput: Must specify electron mass using 'B0'");
  }

  void
  SOL3DElectronTempAtWallCalc::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis->getNumNodes();

    // Figure out what nodes are on the right edge
    int numEdgeQuadNodes = nodalBasis->getNumSurfGaussNodes();
    Lucee::Matrix<double> interpEdgeMatrixLucee(numEdgeQuadNodes, nlocal);
    Lucee::Matrix<double> gaussEdgeOrdinatesLucee(numEdgeQuadNodes, 3);
    // Allocate Eigen matrices
    Eigen::MatrixXd interpEdgeMatrix(numEdgeQuadNodes, nlocal);
    gaussEdgeOrdinates = Eigen::MatrixXd(numEdgeQuadNodes, 3);
    gaussEdgeWeights = std::vector<double>(numEdgeQuadNodes);
    // Get the interpolation matrix for the right edge quadrature points.
    nodalBasis->getSurfUpperGaussQuadData(0, interpEdgeMatrixLucee, gaussEdgeOrdinatesLucee,
      gaussEdgeWeights);

    rightEdgeNodeNums = std::vector<int>(nodalBasis->getNumSurfUpperNodes(0));
    nodalBasis->getSurfUpperNodeNums(0, rightEdgeNodeNums);

    // Copy matrices to eigen objects
    copyLuceeToEigen(interpEdgeMatrixLucee, interpEdgeMatrix);
    copyLuceeToEigen(gaussEdgeOrdinatesLucee, gaussEdgeOrdinates);

    edgeNodeInterpMatrix = Eigen::MatrixXd(numEdgeQuadNodes, rightEdgeNodeNums.size());
    copyLuceeToEigen(interpEdgeMatrixLucee, edgeNodeInterpMatrix);
    
    // Create a lower dimension (2D) interpolation matrix
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
  }

  Lucee::UpdaterStatus
  SOL3DElectronTempAtWallCalc::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Input electron distribution function
    const Lucee::Field<3, double>& distfIn = this->getInp<Lucee::Field<3, double> >(0);
    // Input dynvector of cutoff velocities
    const Lucee::DynVector<double>& cutoffVelocitiesIn = this->getInp<Lucee::DynVector<double> >(1);
    // Output dynvector of electron temperatures at wall
    Lucee::DynVector<double>& elcTempAtWallOut = this->getOut<Lucee::DynVector<double> >(0);

#ifdef HAVE_MPI
// This updater presently will not work in parallel. Eventually, we
// need to fix this, but it requires some MPI-fu to collect the values
// and put them on the appropriate processors. (Ammar Hakim,
// 8/06/2013)
    throw Lucee::Except("SOL3DElectronTempAtWallCalc does not work in parallel!");
#endif

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();

    unsigned nlocal = nodalBasis->getNumNodes();

    double cellCentroid[3];
    int idx[3];
    
    std::vector<double> cutoffVelocities = cutoffVelocitiesIn.getLastInsertedData();

    double rightEdgeCutoffVelocity = cutoffVelocities[1];

    // Need to compute <1>,<vPara>,<mu>,<vPara^2>, where the parallel velocity integration in < >
    // has a lower bound of vCutoff instead of -infinity
    double mom0Total = 0.0;
    double mom1ParaTotal = 0.0;
    double mom1MuTotal = 0.0;
    double mom2ParaTotal = 0.0;

    // start with cell at maximum v_para, looping down to zero until cutoff velocity is reached
    for (int vIndex = globalRgn.getUpper(1)-1; vIndex >= 0; vIndex--)
    { 
      // x coordinate (right edge only)
      idx[0] = globalRgn.getUpper(0)-1;
      // vParallel coordinate
      idx[1] = vIndex;
      // mu coordinate
      idx[2] = globalRgn.getLower(2);
      
      nodalBasis->setIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      double vCellLower = cellCentroid[1] - grid.getDx(1)/2.0;
      double vCellUpper = cellCentroid[1] + grid.getDx(1)/2.0;

      // Check to see if this is the cutoff cell
      if (rightEdgeCutoffVelocity > vCellLower && rightEdgeCutoffVelocity < vCellUpper)
      {
        // Compute cutoff velocity in reference coordinates [-1,1]
        double refCutoff = (rightEdgeCutoffVelocity - cellCentroid[1])/(grid.getDx(1)/2.0);
        // refCoord is also between -1 and +1
        double refCoord[2];
        refCoord[0] = 1;
        std::vector<double> basisAtPoint(nodalBasis->getNumNodes());

        // Make sure we integrate in vParallel from vCutoff to vCellUpper
        for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
        {
          idx[2] = muIndex;
          distfIn.setPtr(distfPtr, idx);
          for (int gaussNodeIndex = 0; gaussNodeIndex < gaussEdgeOrdinates.rows(); gaussNodeIndex++)
          {
            // calculate quadrature points to integrate over
            refCoord[1] = 0.5*(1-refCutoff)*gaussEdgeOrdinates(gaussNodeIndex, 1)
              + 0.5*(refCutoff+1);
            refCoord[2] = gaussEdgeOrdinates(gaussNodeIndex, 2);
            // coordinates in physical space
            double physicalV = cellCentroid[1] + refCoord[1]*grid.getDx(1)/2.0;
            double physicalMu = cellCentroid[2] + refCoord[2]*grid.getDx(2)/2.0;
            // Evaluate all basis functions at coordinate
            nodalBasis->evalBasis(refCoord, basisAtPoint);
            // Evaluate f at this quadrature point
            double fAtPoint = 0.0;
            for (int nodeIndex = 0; nodeIndex < rightEdgeNodeNums.size(); nodeIndex++)
              fAtPoint += distfPtr[rightEdgeNodeNums[nodeIndex]]*basisAtPoint[rightEdgeNodeNums[nodeIndex]];

            // Accumulate result to integrals (gaussEdgeWeights has been premultipled by grid.getDx(1)/2,
            // so we need to multiply by new integration width
            mom0Total += 0.5*(1-refCutoff)*gaussEdgeWeights[gaussNodeIndex]*fAtPoint;
            mom1ParaTotal += 0.5*(1-refCutoff)*gaussEdgeWeights[gaussNodeIndex]*fAtPoint*physicalV;
            mom1MuTotal += 0.5*(1-refCutoff)*gaussEdgeWeights[gaussNodeIndex]*fAtPoint*physicalMu;
            mom2ParaTotal += 0.5*(1-refCutoff)*gaussEdgeWeights[gaussNodeIndex]*fAtPoint*physicalV*physicalV;
          }
        }

        // Get out of for-loop
        break;
      }
      else
      {
        // Integrate over the whole cell
        for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
        {
          idx[2] = muIndex;
          distfIn.setPtr(distfPtr, idx);
          
          // Copy nodes on right edge into a vector
          Eigen::VectorXd rightEdgeData(rightEdgeNodeNums.size());
          
          for (int edgeNodeIndex = 0; edgeNodeIndex < rightEdgeData.rows(); edgeNodeIndex++)
            rightEdgeData(edgeNodeIndex) = distfPtr[rightEdgeNodeNums[edgeNodeIndex]];
          
          // Interpolate nodal data to quadrature points on the edge
          Eigen::VectorXd rightEdgeQuadData = edgeNodeInterpMatrix*rightEdgeData;

          // Compute and accumulate integrals
          for (int quadNodeIndex = 0; quadNodeIndex < rightEdgeQuadData.rows(); quadNodeIndex++)
          {
            double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
            double physicalMu = cellCentroid[2] + gaussEdgeOrdinates(quadNodeIndex,2)*grid.getDx(2)/2.0;

            mom0Total += gaussEdgeWeights[quadNodeIndex]*rightEdgeQuadData(quadNodeIndex);
            mom1ParaTotal += gaussEdgeWeights[quadNodeIndex]*physicalV*rightEdgeQuadData(quadNodeIndex);
            mom1MuTotal += gaussEdgeWeights[quadNodeIndex]*physicalMu*rightEdgeQuadData(quadNodeIndex);
            mom2ParaTotal += gaussEdgeWeights[quadNodeIndex]*physicalV*physicalV*rightEdgeQuadData(quadNodeIndex);
          }
        }
      }
    }

    double uParallel = mom1ParaTotal/mom0Total;
    std::vector<double> data(1);
    // compute electron temperature in units of eV
    data[0] = 2.0/3.0*(0.5*elcMass*mom2ParaTotal/mom0Total
        + B0*mom1MuTotal/mom0Total -
        0.5*elcMass*rightEdgeCutoffVelocity*rightEdgeCutoffVelocity)/eV;
    
    // Put data into the DynVector
    elcTempAtWallOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  SOL3DElectronTempAtWallCalc::declareTypes()
  {
    // Input: Full 3d distribution function
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: DynVector containing cutoff velocities
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // Output: DynVector containing electron temperature at wall vs. time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  SOL3DElectronTempAtWallCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

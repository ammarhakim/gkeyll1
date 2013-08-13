/**
 * @file	LcDistFuncReflectionBcUpdater.cpp
 *
 * @brief	Applies particle refection BCs to distribution function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncReflectionBcUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;
  static const unsigned LC_BOTH_EDGES = 2;

  const char *DistFuncReflectionBcUpdater::id = "DistFuncReflectionBc";

  DistFuncReflectionBcUpdater::DistFuncReflectionBcUpdater()
  {
  }

  void
  DistFuncReflectionBcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("DistFuncReflectionBcUpdater::readInput: Must specify element to use using 'basis'");

    applyLeftEdge = applyRightEdge = false;
    std::string edgeStr = tbl.getString("edge");
    if (edgeStr == "lower")
      applyLeftEdge = true;
    else if (edgeStr == "upper")
      applyRightEdge = true;
    else if (edgeStr == "both")
      applyLeftEdge = applyRightEdge = true;
  }

  void
  DistFuncReflectionBcUpdater::initialize()
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

    // Tolerance for finding cutoff velocities
    cutoffTolerance = 0.0001;
  }

  Lucee::UpdaterStatus
  DistFuncReflectionBcUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Momentum flux for ions
    const Lucee::Field<1, double>& mom1 = this->getInp<Lucee::Field<1, double> >(0);
    // Output distribution function
    Lucee::Field<2, double>& distf = this->getOut<Lucee::Field<2, double> >(0);
    // Output cutoff velocities vs time
    Lucee::DynVector<double>& cutoffVelocities = this->getOut<Lucee::DynVector<double> >(1);

#ifdef HAVE_MPI
// This updater presently will not work in parallel. Eventually, we
// need to fix this, but it requires some MPI-fu to collect the values
// and put them on the appropriate processors. (Ammar Hakim,
// 8/06/2013)
    throw Lucee::Except("DistFuncReflectionBcUpdater does not work in parallel!");
#endif

    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();
    Lucee::ConstFieldPtr<double> mom1Ptr = mom1.createConstPtr();
    Lucee::FieldPtr<double> sknPtr = distf.createPtr(); // for skin-cell
    Lucee::FieldPtr<double> gstPtr = distf.createPtr(); // for ghost-cell

    unsigned numNodes = nodalBasis->getNumNodes();

    std::vector<double> data(2);
    // Initialize data to 0's. Convention: 0 = left, 1 = right
    data[0] = 0.0;
    data[1] = 0.0;

    if (applyRightEdge)
    {
      int ix = globalRgn.getUpper(0)-1; // right skin cell x index
      double totalFluxAlongEdge = 0.0; 
      double cellCentroid[3];
      int idx[2];
      bool foundCutoffVelocity = false;

      // Find flux at right-most node
      mom1.setPtr(mom1Ptr, ix);
      // Assume that the location of the right-most node in 1-D
      // is the number of nodes along an edge in 2-D minus 1
      double totalIonFlux = mom1Ptr[rightEdgeNodeNums.size()-1];

      for (int js=globalRgn.getUpper(1)-1, jg=0; js>=0; --js, ++jg)
      {
        nodalBasis->setIndex(ix, js);
      
        idx[0] = ix;
        idx[1] = js;
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't need to loop over left-moving cells
        if (cellCentroid[1] < 0.0)
          break;

        distf.setPtr(sknPtr, ix, js);
        distf.setPtr(gstPtr, ix+1, jg);

        // Copy nodes on right edge into a vector
        Eigen::VectorXd rightEdgeData(rightEdgeNodeNums.size());
        for (int edgeNodeIndex = 0; edgeNodeIndex < rightEdgeData.rows(); edgeNodeIndex++)
          rightEdgeData(edgeNodeIndex) = sknPtr[rightEdgeNodeNums[edgeNodeIndex]];
        
        // Interpolate nodal data to quadrature points on the edge
        Eigen::VectorXd rightEdgeQuadData = edgeNodeInterpMatrix*rightEdgeData;
        
        // Integrate v*f over entire cell using gaussian quadrature
        double fluxInEntireCell = 0.0;
        for (int quadNodeIndex = 0; quadNodeIndex < rightEdgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
          fluxInEntireCell += gaussEdgeWeights[quadNodeIndex]*physicalV*rightEdgeQuadData(quadNodeIndex);
        }

        if (foundCutoffVelocity == false)
        {
          if (totalFluxAlongEdge + fluxInEntireCell < totalIonFlux)
          {
            totalFluxAlongEdge += fluxInEntireCell;

            // Set to no inflow condition then go on to next cell
            for (int componentIndex = 0; componentIndex < numNodes; componentIndex++)
              gstPtr[componentIndex] = 0.0;
          }
          else
          {
            // This is the cell that contains the cutoff velocity
            
            // Figure out fraction of cell that contains the excess flux
            // (Flux over what is needed for equivalence with Gamma_i)
            double excessFraction = (fluxInEntireCell + totalFluxAlongEdge - totalIonFlux)/fluxInEntireCell;

            // cutoffGuess is between -deltaV/2 and +deltaV/2
            double cutoffGuess = - grid.getDx(1)/2.0 + excessFraction*grid.getDx(1);
            double exactResult = totalIonFlux - totalFluxAlongEdge;
            double cellWidth = grid.getDx(1)/2.0;
  
            if (exactResult != 0)
              cutoffGuess = findRightCutoffVelocity(sknPtr, cutoffGuess, exactResult, cellWidth, cellCentroid, cutoffTolerance);

            foundCutoffVelocity = true;
            data[1] = cellCentroid[1] + cutoffGuess;
            
            // Scale all values by appropriate fraction to conserve flux
            // Copy data into ghost after rotating skin cell data by 180 degrees
            for (int k = 0; k < numNodes; k++)
              gstPtr[k] = excessFraction*sknPtr[rotMap[k]];
          }
        }
        else
        {
          // We have already iterated past and found the cutoff velocity.
          // Copy data into ghost after rotating skin cell data by 180 degrees
          for (int k = 0; k < numNodes; k++)
            gstPtr[k] = sknPtr[rotMap[k]];
        }
      }
    }
    
    if (applyLeftEdge)
    {
      int ix = globalRgn.getLower(0); // left skin cell x index
      double totalFluxAlongEdge = 0.0; 
      double cellCentroid[3];
      int idx[2];
      bool foundCutoffVelocity = false;

      // Find flux at left-most node
      mom1.setPtr(mom1Ptr, ix);
      // Assume that the location of the left-most node in 1-D is 0
      // This should be a negative number
      double totalIonFlux = mom1Ptr[0];

      for (int js=0, jg=globalRgn.getUpper(1)-1; jg>=0; ++js, --jg)
      {
        nodalBasis->setIndex(ix, js);

        distf.setPtr(sknPtr, ix, js);
        distf.setPtr(gstPtr, ix-1, jg);
      
        idx[0] = ix;
        idx[1] = js;
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't need to loop over right-moving cells
        if (cellCentroid[1] > 0.0)
          break;

        // Copy nodes on left edge into a vector
        Eigen::VectorXd leftEdgeData(leftEdgeNodeNums.size());
        for (int edgeNodeIndex = 0; edgeNodeIndex < leftEdgeData.rows(); edgeNodeIndex++)
          leftEdgeData(edgeNodeIndex) = sknPtr[leftEdgeNodeNums[edgeNodeIndex]];
        
        // Interpolate nodal data to quadrature points on the edge
        Eigen::VectorXd leftEdgeQuadData = edgeNodeInterpMatrix*leftEdgeData;
        
        // Integrate v*f over entire cell using gaussian quadrature
        double fluxInEntireCell = 0.0;
        for (int quadNodeIndex = 0; quadNodeIndex < leftEdgeQuadData.rows(); quadNodeIndex++)
        {
          double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
          fluxInEntireCell += gaussEdgeWeights[quadNodeIndex]*physicalV*leftEdgeQuadData(quadNodeIndex);
        }

        if (foundCutoffVelocity == false)
        {
          if (totalFluxAlongEdge + fluxInEntireCell > totalIonFlux)
          {
            totalFluxAlongEdge += fluxInEntireCell;

            // Set to no inflow condition then go on to next cell
            for (int componentIndex = 0; componentIndex < numNodes; componentIndex++)
              gstPtr[componentIndex] = 0.0;
          }
          else
          {
            // This is the cell that contains the cutoff velocity
            
            // Figure out fraction of cell that contains the excess flux
            // (Flux over what is needed for equivalence with Gamma_i)
            double excessFraction = (fluxInEntireCell + totalFluxAlongEdge - totalIonFlux)/fluxInEntireCell;

            // cutoffGuess is between -deltaV/2 and +deltaV/2
            double cutoffGuess = - grid.getDx(1)/2.0 + excessFraction*grid.getDx(1);
            double exactResult = totalIonFlux - totalFluxAlongEdge;
            double cellWidth = grid.getDx(1)/2.0;

            if (exactResult != 0)
              cutoffGuess = findLeftCutoffVelocity(sknPtr, cutoffGuess, exactResult, cellWidth, cellCentroid, cutoffTolerance);
            
            foundCutoffVelocity = true;
            data[0] = cellCentroid[1] + cutoffGuess;
            
            // Scale all values by appropriate fraction to conserve flux
            // Copy data into ghost after rotating skin cell data by 180 degrees
            for (int k = 0; k < numNodes; k++)
              gstPtr[k] = excessFraction*sknPtr[rotMap[k]];
          }
        }
        else
        {
          // We have already iterated past and found the cutoff velocity.
          // Copy data into ghost after rotating skin cell data by 180 degrees
          for (int k = 0; k < numNodes; k++)
            gstPtr[k] = sknPtr[rotMap[k]];
        }
      }
    }

    // Put data into the DynVector
    cutoffVelocities.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncReflectionBcUpdater::declareTypes()
  {
    // Momentum flux for ion
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Distribution function output
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    // Output cutoff velocities
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  DistFuncReflectionBcUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

  double
  DistFuncReflectionBcUpdater::findRightCutoffVelocity(const Lucee::ConstFieldPtr<double>& searchFld,
    const double initialGuess, const double exactResult, const double cellWidth, const double* cellCentroid, const double tol)
  {
    double refCoord[2];
    refCoord[0] = 1;

    std::vector<double> basisAtPoint(nodalBasis->getNumNodes());

    double b = cellWidth;
    double cutoffGuess;
    double nextCutoffGuess = initialGuess;
    
    double upperBound = cellWidth;
    double lowerBound = -cellWidth;
    
    double relError; 

    do
    {
      cutoffGuess = nextCutoffGuess;

      double integralResult = 0.0;
      for (int gaussNodeIndex = 0; gaussNodeIndex < gaussEdgeOrdinates.rows(); gaussNodeIndex++)
      {
        // physicalCoord is between -deltaV/2 and +deltaV/2
        double physicalCoord = 0.5*(b-cutoffGuess)*gaussEdgeOrdinates(gaussNodeIndex, 1)
          + 0.5*(b + cutoffGuess);
        refCoord[1] = physicalCoord/cellWidth;
        
        nodalBasis->evalBasis(refCoord, basisAtPoint);
        // Evaluate f at this location
        double fAtPoint = 0.0;
        // Loop over 1-D basis functions
        for (int nodeIndex = 0; nodeIndex < rightEdgeNodeNums.size(); nodeIndex++)
          fAtPoint += searchFld[rightEdgeNodeNums[nodeIndex]]*basisAtPoint[rightEdgeNodeNums[nodeIndex]];

        integralResult += 0.5*(b - cutoffGuess)*gaussEdgeWeights[gaussNodeIndex]/cellWidth*(cellCentroid[1] + physicalCoord)*fAtPoint;
      }
      
      relError = (integralResult - exactResult)/exactResult;
      
      if (relError < 0)
        upperBound = cutoffGuess;
      else
        lowerBound = cutoffGuess;

      nextCutoffGuess = 0.5*(lowerBound + upperBound);
    }
    while (fabs(relError) > tol);

    return cutoffGuess;
  }

  double
  DistFuncReflectionBcUpdater::findLeftCutoffVelocity(const Lucee::ConstFieldPtr<double>& searchFld,
    const double initialGuess, const double exactResult, const double cellWidth, const double* cellCentroid, const double tol)
  {
    double refCoord[2];
    refCoord[0] = -1;

    std::vector<double> basisAtPoint(nodalBasis->getNumNodes());

    double b = -cellWidth;
    double cutoffGuess;
    double nextCutoffGuess = initialGuess;
    
    double upperBound = cellWidth;
    double lowerBound = -cellWidth;
    
    double relError; 

    do
    {
      cutoffGuess = nextCutoffGuess;

      double integralResult = 0.0;
      for (int gaussNodeIndex = 0; gaussNodeIndex < gaussEdgeOrdinates.rows(); gaussNodeIndex++)
      {
        // physicalCoord is between -deltaV/2 and +deltaV/2
        double physicalCoord = 0.5*(cutoffGuess-b)*gaussEdgeOrdinates(gaussNodeIndex, 1)
          + 0.5*(cutoffGuess + b);
        refCoord[1] = physicalCoord/cellWidth;
        
        nodalBasis->evalBasis(refCoord, basisAtPoint);
        // Evaluate f at this location
        double fAtPoint = 0.0;
        // Loop over 1-D basis functions
        for (int nodeIndex = 0; nodeIndex < leftEdgeNodeNums.size(); nodeIndex++)
          fAtPoint += searchFld[leftEdgeNodeNums[nodeIndex]]*basisAtPoint[leftEdgeNodeNums[nodeIndex]];

        integralResult += 0.5*(cutoffGuess - b)*gaussEdgeWeights[gaussNodeIndex]/cellWidth*(cellCentroid[1] + physicalCoord)*fAtPoint;
      }
      
      relError = (integralResult - exactResult)/exactResult;
      
      if (relError > 0)
        upperBound = cutoffGuess;
      else
        lowerBound = cutoffGuess;

      nextCutoffGuess = 0.5*(lowerBound + upperBound);
    }
    while (fabs(relError) > tol);

    return cutoffGuess;
  }
}

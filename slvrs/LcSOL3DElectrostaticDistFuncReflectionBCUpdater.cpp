/**
 * @file	LcDistFuncReflectionBcUpdater.cpp
 *
 * @brief	Applies particle refection BCs to (electron) distribution function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOL3DElectrostaticDistFuncReflectionBCUpdater.h>
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

  const char *SOL3DElectrostaticDistFuncReflectionBCUpdater::id = "SOL3DElectrostaticDistFuncReflectionBc";

  SOL3DElectrostaticDistFuncReflectionBCUpdater::SOL3DElectrostaticDistFuncReflectionBCUpdater()
  {
  }

  void
  SOL3DElectrostaticDistFuncReflectionBCUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOL3DElectrostaticDistFuncReflectionBCUpdater::readInput: Must specify element to use using 'basis'");

    // get hold of 1D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except(
        "SOL3DElectrostaticDistFuncReflectionBCUpdater::readInput: Must specify 1D element to use using 'basis1d'");

    applyLeftEdge = applyRightEdge = false;
    std::string edgeStr = tbl.getString("edge");
    if (edgeStr == "lower")
      applyLeftEdge = true;
    else if (edgeStr == "upper")
      applyRightEdge = true;
    else if (edgeStr == "both")
      applyLeftEdge = applyRightEdge = true;

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    if (tbl.hasBool("computeCutoffVelocities"))
      computeCutoffVelocities = tbl.getBool("computeCutoffVelocities");
    else computeCutoffVelocities = true;
  }

  void
  SOL3DElectrostaticDistFuncReflectionBCUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<unsigned> yRef(nlocal), xRef(nlocal);

    // reflection mapping
    rotMap.resize(nlocal);
    nodalBasis->getUpperReflectingBcMapping(0, yRef);
    nodalBasis->getLowerReflectingBcMapping(1, xRef);
    for (int i = 0; i < nlocal; i++)
      rotMap[i] = xRef[yRef[i]];


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
    leftEdgeNodeNums = std::vector<int>(nodalBasis->getNumSurfLowerNodes(0));
    nodalBasis->getSurfUpperNodeNums(0, rightEdgeNodeNums);
    nodalBasis->getSurfLowerNodeNums(0, leftEdgeNodeNums);

    // Copy matrices to eigen objects
    copyLuceeToEigen(interpEdgeMatrixLucee, interpEdgeMatrix);
    copyLuceeToEigen(gaussEdgeOrdinatesLucee, gaussEdgeOrdinates);
    edgeNodeInterpMatrix = Eigen::MatrixXd(numEdgeQuadNodes, rightEdgeNodeNums.size());
    
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

    // Temporary testing code
    gaussEdgeWeightsLeftEdge = std::vector<double>(numEdgeQuadNodes);
    gaussEdgeOrdinatesLeftEdge = Eigen::MatrixXd(numEdgeQuadNodes, 3);
    nodalBasis->getSurfLowerGaussQuadData(0, interpEdgeMatrixLucee, gaussEdgeOrdinatesLucee,
      gaussEdgeWeightsLeftEdge);
    copyLuceeToEigen(interpEdgeMatrixLucee, interpEdgeMatrix);
    copyLuceeToEigen(gaussEdgeOrdinatesLucee, gaussEdgeOrdinatesLeftEdge);
    edgeNodeInterpMatrixLeftEdge = Eigen::MatrixXd(numEdgeQuadNodes, leftEdgeNodeNums.size());
    // Create a lower dimension (2D) interpolation matrix
    for (int nodeIndex = 0; nodeIndex < numEdgeQuadNodes; nodeIndex++)
    {
      // At each quadrature node, copy basis function evaluations for
      // those basis functions associated with the nodes on the (right) edge
      for (int basisIndex = 0; basisIndex < leftEdgeNodeNums.size(); basisIndex++)
      {
        edgeNodeInterpMatrixLeftEdge(nodeIndex, basisIndex) = interpEdgeMatrix(nodeIndex, 
          leftEdgeNodeNums[basisIndex]);
      }
    }

    // Tolerance for finding cutoff velocities
    cutoffTolerance = 1.0e-4;
  }

  Lucee::UpdaterStatus
  SOL3DElectrostaticDistFuncReflectionBCUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::DynVector<double>& momentsAtEdgesIonIn = this->getInp<Lucee::DynVector<double> >(0);
    // Output distribution function
    Lucee::Field<3, double>& distf = this->getOut<Lucee::Field<3, double> >(0);
    // Output cutoff velocities vs time
    Lucee::DynVector<double>& cutoffVelocities = this->getOut<Lucee::DynVector<double> >(1);

#ifdef HAVE_MPI
// This updater presently will not work in parallel. Eventually, we
// need to fix this, but it requires some MPI-fu to collect the values
// and put them on the appropriate processors. (Ammar Hakim,
// 8/06/2013)
    throw Lucee::Except("SOL3DElectrostaticDistFuncReflectionBCUpdater does not work in parallel!");
#endif

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    const std::vector<double> momentsAtEdgesIon = momentsAtEdgesIonIn.getLastInsertedData();
    Lucee::FieldPtr<double> sknPtr = distf.createPtr(); // for skin-cell
    Lucee::FieldPtr<double> gstPtr = distf.createPtr(); // for ghost-cell

    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned nlocal1d = nodalBasis1d->getNumNodes();

    std::vector<double> data(2);
    // Initialize data to 0's. Convention: 0 = left edge, 1 = right edge
    data[0] = 0.0;
    data[1] = 0.0;

    if (applyRightEdge)
    {
      int ix = globalRgn.getUpper(0)-1; // right skin cell x index
      double cellCentroid[3];
      int idx[3];
      bool foundCutoffVelocity = false;
      double totalFluxAlongEdge = 0.0; 

      // Get ion flux at right edge
      double totalIonFlux = momentsAtEdgesIon[0];

      // start with cell at maximum v_para, looping down to zero until desired flux is reached
      for (int js=globalRgn.getUpper(1)-1, jg=0; js>=0; --js, ++jg)
      {
        idx[0] = ix;
        idx[1] = js;
        idx[2] = globalRgn.getLower(2);
        
        nodalBasis->setIndex(idx);
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);
        
        // Don't need to loop over left-moving cells
        if (cellCentroid[1] < 0.0)
          break;
        
        double fluxInEntireCell = 0.0;
        // At this value of v_parallel, need to integrate flux over all mu
        for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
        {
          idx[2] = muIndex;
          distf.setPtr(sknPtr, ix, js, muIndex);

          // Copy nodes on right edge into a vector
          Eigen::VectorXd rightEdgeData(rightEdgeNodeNums.size());
          
          for (int edgeNodeIndex = 0; edgeNodeIndex < rightEdgeData.rows(); edgeNodeIndex++)
            rightEdgeData(edgeNodeIndex) = sknPtr[rightEdgeNodeNums[edgeNodeIndex]];
          
          // Interpolate nodal data to quadrature points on the edge
          Eigen::VectorXd rightEdgeQuadData = edgeNodeInterpMatrix*rightEdgeData;
        
          // Integrate v*f over entire cell using gaussian quadrature
          for (int quadNodeIndex = 0; quadNodeIndex < rightEdgeQuadData.rows(); quadNodeIndex++)
          {
            double physicalV = cellCentroid[1] + gaussEdgeOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
            fluxInEntireCell += scaleFactor*gaussEdgeWeights[quadNodeIndex]*physicalV*
              rightEdgeQuadData(quadNodeIndex);
          }
        }

        if (foundCutoffVelocity == false)
        {
          if (totalFluxAlongEdge + fluxInEntireCell < totalIonFlux)
          {
            // Accumulate to total since cutoff velocity is not in this cell
            totalFluxAlongEdge += fluxInEntireCell;
            // Set to no inflow condition then go on to next cell
            for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
            {
              distf.setPtr(gstPtr, ix+1, jg, muIndex);
              
              for (int k = 0; k < nlocal; k++)
                gstPtr[k] = 0.0;
            }
          }
          else
          {
            // This is the cell that contains the cutoff velocity
            
            // Figure out fraction of cell that contains the excess flux
            // (Flux over what is needed for equivalence with Gamma_i)
            double excessFraction = (fluxInEntireCell + totalFluxAlongEdge - totalIonFlux)/fluxInEntireCell;
            //std::cout << "excessFraction (L) = " << excessFraction << std::endl;
            // Search for cutoff velocity if needed
            if (computeCutoffVelocities == true)
            {
              // desired result for integrating from flux vCutoff to vUpper in this cell
              double exactResult = totalIonFlux - totalFluxAlongEdge;
              double cellWidth = grid.getDx(1)/2.0;
              double cutoffGuess;

              if (exactResult != 0.0 && totalIonFlux != 0.0)
              {
                double refCoord[3];
                refCoord[0] = 1;

                std::vector<double> basisAtPoint(rightEdgeNodeNums.size());

                // upper boundary for integration
                double b = grid.getDx(1)/2.0;
                // cutoffGuess is between -dV/2 and +dV/2
                double nextCutoffGuess = -grid.getDx(1)/2.0 + excessFraction*grid.getDx(1);

                //double nextCutoffGuess = -grid.getDx(1)/fluxInEntireCell*( totalIonFlux - totalFluxAlongEdge - fluxInEntireCell/2.0 );
                
                double upperBound = grid.getDx(1)/2.0;
                double lowerBound = -grid.getDx(1)/2.0;
                
                double relError; 

                int iterCount = 0;

                do
                {
                  cutoffGuess = nextCutoffGuess;

                  if (upperBound == lowerBound || iterCount > 100)
                    break;

                  double integralResult = 0.0;
                  // Integrate flux from vCutoff to vUpper
                  for (int gaussNodeIndex = 0; gaussNodeIndex < gaussEdgeOrdinates.rows(); gaussNodeIndex++)
                  {
                    // physicalCoord (vParallel) is between -deltaV/2 and +deltaV/2
                    double physicalCoord = 0.5*(b-cutoffGuess)*gaussEdgeOrdinates(gaussNodeIndex, 1)
                      + 0.5*(b + cutoffGuess);
                    // refCoord is between -1 and +1
                    refCoord[1] = physicalCoord/(grid.getDx(1)/2.0);
                    refCoord[2] = gaussEdgeOrdinates(gaussNodeIndex, 2);
                    // Evaluate all basis functions at coordinate
                    nodalBasis->evalBasis(refCoord, basisAtPoint, rightEdgeNodeNums);
                    
                    for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
                    {
                      distf.setPtr(sknPtr, ix, js, muIndex);
                      // Evaluate f at this quadrature point
                      double fAtPoint = 0.0;
                      // Loop over 2D basis functions
                      for (int nodeIndex = 0; nodeIndex < rightEdgeNodeNums.size(); nodeIndex++)
                        fAtPoint += sknPtr[rightEdgeNodeNums[nodeIndex]]*basisAtPoint[nodeIndex];

                      // Need to divide by gaussEdgeWeights by cell width in V_parallel then multiply by
                      // integration width
                      integralResult += scaleFactor*(b - cutoffGuess)*gaussEdgeWeights[gaussNodeIndex]/grid.getDx(1)*fAtPoint*(cellCentroid[1] + physicalCoord);
                    }
                  }
                  
                  relError = (integralResult - exactResult)/exactResult;

                  if (iterCount > 50)
                  {
                    std::cout << "iterCount = " << iterCount << std::endl;
                    std::cout << "relError = " << fabs(relError) << std::endl;
                    std::cout << "integralResult = " << integralResult << std::endl;
                    std::cout << "exactResult = " << exactResult << std::endl;
                    std::cout << "upperBound = " << upperBound/cellWidth << std::endl;
                    std::cout << "lowerBound = " << lowerBound/cellWidth << std::endl;
                  }

                  if (relError < 0)
                    upperBound = cutoffGuess;
                  else
                    lowerBound = cutoffGuess;

                  nextCutoffGuess = 0.5*(lowerBound + upperBound);
                  iterCount++;
                }
                while (fabs(relError) > cutoffTolerance);
                //std::cout << "Iterations (R) = " << iterCount << std::endl;
              }

              // store cutoff velocity
              data[1] = cellCentroid[1] + cutoffGuess;
            }
            
            foundCutoffVelocity = true;
            // Scale all values by appropriate fraction to conserve flux
            // Copy data into ghost after rotating skin cell data by 180 degrees
            for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
            {
              distf.setPtr(sknPtr, ix, js, muIndex);
              distf.setPtr(gstPtr, ix+1, jg, muIndex);

              for (int k = 0; k < nlocal; k++)
                gstPtr[k] = excessFraction*sknPtr[rotMap[k]];
            }
          }
        }
        else
        {
          // We have already iterated past and found the cutoff velocity.
          // Copy data into ghost after rotating skin cell data by 180 degrees
          for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
          {
            distf.setPtr(gstPtr, ix+1, jg, muIndex);
            distf.setPtr(sknPtr, ix, js, muIndex);

            for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
              gstPtr[componentIndex] = sknPtr[rotMap[componentIndex]];
          }
        }
      }
    }
    
    
    if (applyLeftEdge)
    {
      int ix = globalRgn.getLower(0); // left skin cell x index
      double cellCentroid[3];
      int idx[3];
      bool foundCutoffVelocity = false;
      double totalFluxAlongEdge = 0.0; 

      // Find flux at left edge
      double totalIonFlux = momentsAtEdgesIon[3];

      for (int js=0, jg=globalRgn.getUpper(1)-1; jg>=0; ++js, --jg)
      {
        idx[0] = ix;
        idx[1] = js;
        idx[2] = globalRgn.getLower(2);

        nodalBasis->setIndex(idx);
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Don't need to loop over right-moving cells
        if (cellCentroid[1] > 0.0)
          break;

        double fluxInEntireCell = 0.0;
        // At this value of v_parallel, need to integrate flux over all mu
        for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
        {
          idx[2] = muIndex;
          distf.setPtr(sknPtr, ix, js, muIndex);
          
          // Copy nodes on left edge into a vector
          Eigen::VectorXd leftEdgeData(leftEdgeNodeNums.size());
          
          for (int edgeNodeIndex = 0; edgeNodeIndex < leftEdgeData.rows(); edgeNodeIndex++)
            leftEdgeData(edgeNodeIndex) = sknPtr[leftEdgeNodeNums[edgeNodeIndex]]; 
          
          // Interpolate nodal data to quadrature points on the edge
          Eigen::VectorXd leftEdgeQuadData = edgeNodeInterpMatrixLeftEdge*leftEdgeData;
        
          // Integrate v*f over entire cell using gaussian quadrature
          for (int quadNodeIndex = 0; quadNodeIndex < leftEdgeQuadData.rows(); quadNodeIndex++)
          {
            double physicalV = cellCentroid[1] + gaussEdgeOrdinatesLeftEdge(quadNodeIndex,1)*grid.getDx(1)/2.0;
            fluxInEntireCell += scaleFactor*gaussEdgeWeightsLeftEdge[quadNodeIndex]*physicalV*
              leftEdgeQuadData(quadNodeIndex);
          }
        }

        if (foundCutoffVelocity == false)
        {
          if (totalFluxAlongEdge + fluxInEntireCell > totalIonFlux)
          {
            totalFluxAlongEdge += fluxInEntireCell;

            // Set to no inflow condition then go on to next cell
            for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
            {
              distf.setPtr(gstPtr, ix-1, jg, muIndex);
              
              for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
                gstPtr[componentIndex] = 0.0;
            }
          }
          else
          {
            // This is the cell that contains the cutoff velocity
            
            // Figure out fraction of cell that contains the excess flux
            // (Flux over what is needed for equivalence with Gamma_i)
            double excessFraction = (fluxInEntireCell + totalFluxAlongEdge - totalIonFlux)/fluxInEntireCell;
            //std::cout << "excessFraction (R) = " << excessFraction << std::endl;

            // Search for cutoff velocity if needed
            if (computeCutoffVelocities == true)
            {
              // desired result for integrating from flux vCutoff to vUpper in this cell
              double exactResult = totalIonFlux - totalFluxAlongEdge;
              double cellWidth = grid.getDx(1)/2.0;
              double cutoffGuess;

              if (exactResult != 0.0 && totalIonFlux != 0.0)
              {
                double refCoord[3];
                refCoord[0] = -1;

                std::vector<double> basisAtPoint(leftEdgeNodeNums.size());

                // upper boundary for integration
                double b = -grid.getDx(1)/2.0;
                // cutoffGuess is between -dV/2 and +dV/2
                double nextCutoffGuess = grid.getDx(1)/2.0 - excessFraction*grid.getDx(1);
                // Initial guess will be negative of cutoff velocity found on right edge
                if(applyRightEdge == true)
                  nextCutoffGuess = -data[1] - cellCentroid[1];
                
                double upperBound = grid.getDx(1)/2.0;
                double lowerBound = -grid.getDx(1)/2.0;
                
                double relError; 

                int iterCount = 0;

                do
                {
                  cutoffGuess = nextCutoffGuess;

                  if (upperBound == lowerBound || iterCount > 100)
                    break;

                  double integralResult = 0.0;
                  // Integrate flux from vCutoff to vUpper
                  for (int gaussNodeIndex = 0; gaussNodeIndex < gaussEdgeOrdinates.rows(); gaussNodeIndex++)
                  {
                    // physicalCoord is between -deltaV/2 and +deltaV/2
                    double physicalCoord = 0.5*(cutoffGuess-b)*gaussEdgeOrdinates(gaussNodeIndex, 1)
                      + 0.5*(cutoffGuess + b);
                    // refCoord is between -1 and +1
                    refCoord[1] = physicalCoord/cellWidth;
                    refCoord[2] = gaussEdgeOrdinates(gaussNodeIndex, 2);
                    
                    // Evaluate all basis functions at coordinate
                    nodalBasis->evalBasis(refCoord, basisAtPoint, leftEdgeNodeNums);
                    for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
                    {
                      distf.setPtr(sknPtr, ix, js, muIndex);
                      // Evaluate f at this quadrature point
                      double fAtPoint = 0.0;
                      // Loop over 2D basis functions
                      for (int nodeIndex = 0; nodeIndex < leftEdgeNodeNums.size(); nodeIndex++)
                        fAtPoint += sknPtr[leftEdgeNodeNums[nodeIndex]]*basisAtPoint[nodeIndex];

                      // Need to divide by gaussEdgeWeights by cell width in V_parallel then multiply by
                      // integration width
                      integralResult += scaleFactor*(cutoffGuess-b)*gaussEdgeWeights[gaussNodeIndex]/grid.getDx(1)*fAtPoint*(cellCentroid[1] + physicalCoord);
                    }
                  }
                  
                  relError = (integralResult - exactResult)/exactResult;

                  if (iterCount > 50)
                  {
                    std::cout << "iterCount = " << iterCount << std::endl;
                    std::cout << "relError = " << fabs(relError) << std::endl;
                    std::cout << "integralResult = " << integralResult << std::endl;
                    std::cout << "exactResult = " << exactResult << std::endl;
                    std::cout << "upperBound = " << upperBound/cellWidth << std::endl;
                    std::cout << "lowerBound = " << lowerBound/cellWidth << std::endl;
                  }

                  if (relError > 0)
                    upperBound = cutoffGuess;
                  else
                    lowerBound = cutoffGuess;

                  nextCutoffGuess = 0.5*(lowerBound + upperBound);
                  iterCount++;
                }
                while (fabs(relError) > cutoffTolerance);
                //std::cout << "Iterations (L) = " << iterCount << std::endl;
              }

              // store cutoff velocity
              data[0] = cellCentroid[1] + cutoffGuess;
            }

            foundCutoffVelocity = true;
            
            // Scale all values by appropriate fraction to conserve flux
            // Copy data into ghost after rotating skin cell data by 180 degrees
            for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
            {
              distf.setPtr(sknPtr, ix, js, muIndex);
              distf.setPtr(gstPtr, ix-1, jg, muIndex);

              for (int k = 0; k < nlocal; k++)
                gstPtr[k] = excessFraction*sknPtr[rotMap[k]];
            }
          }
        }
        else
        {
          // We have already iterated past and found the cutoff velocity.
          // Copy data into ghost after rotating skin cell data by 180 degrees
          // Set to no inflow condition then go on to next cell
          for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
          {
            distf.setPtr(sknPtr, ix, js, muIndex);
            distf.setPtr(gstPtr, ix-1, jg, muIndex);

            for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
              gstPtr[componentIndex] = sknPtr[rotMap[componentIndex]];
          }
        }
      }
    }

    // Put data into the DynVector
    if (computeCutoffVelocities == true)
      cutoffVelocities.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  SOL3DElectrostaticDistFuncReflectionBCUpdater::declareTypes()
  {
    // Moments of the ion distribution function on both edges
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // Distribution function output (electrons)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // Output cutoff velocities
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  SOL3DElectrostaticDistFuncReflectionBCUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

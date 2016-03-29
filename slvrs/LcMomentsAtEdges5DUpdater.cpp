/**
 * @file	LcMomentsAtEdges5DUpdater.cpp
 *
 * @brief	Computes several parallel velocity moments of the distribution function at both edges.
 * Used for 5D SOL problem
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMomentsAtEdges5DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

//etc includes
#include <quadrule.hpp>

namespace Lucee
{
  const char *MomentsAtEdges5DUpdater::id = "MomentsAtEdges5D";

  MomentsAtEdges5DUpdater::MomentsAtEdges5DUpdater()
  {
  }

  void
  MomentsAtEdges5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("MomentsAtEdges5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("MomentsAtEdges5DUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("MomentsAtEdges5DUpdater::readInput: Must specify basis function order using 'polyOrder'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    // Indicates if we will integrate in ghost cells in position space
    integrateGhosts = false;
    if (tbl.hasBool("integrateGhosts"))
      integrateGhosts = tbl.getBool("integrateGhosts");
  }

  void
  MomentsAtEdges5DUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis5d->getNumNodes();
    std::vector<unsigned> zRef(nlocal), vRef(nlocal);

    // Get a copy of the nodal coordinates
    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 5);
    nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
    Eigen::MatrixXd nodeCoords(nlocal, 5);
    copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Used to figure out which nodes share the same location in configuration space
    double dxMin = grid.getDx(0);
    for (int d = 1; d < 3; d++)
      dxMin = std::min(dxMin, grid.getDx(d));

    // Find all nodes that share the same location as node zero in configuration space
    nodalStencil = std::vector<int>(nlocal);
    int stencilIndex = 0;
    for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
    {
      if (sameConfigCoords(0, nodeIndex, dxMin, nodeCoords) == true)
      {
        nodalStencil[stencilIndex] = nodeIndex;
        stencilIndex++;
      }
    }
    nodalStencil.resize(stencilIndex);

    // Construct a quadrature matrix to integrate over entire (v,mu) cell
    int integrationDegree = polyOrder + 1; // (need to integrate basis function times one degree in v)
    unsigned numGaussPoints1d = (unsigned)((integrationDegree+1)/2.0 + 0.5);
    std::vector<double> gaussPoints1d(numGaussPoints1d);
    std::vector<double> gaussWeights1d(numGaussPoints1d);
    legendre_set(numGaussPoints1d, &gaussPoints1d[0], &gaussWeights1d[0]);

    int totalSurfQuadPoints = numGaussPoints1d*numGaussPoints1d;
    Eigen::MatrixXd gaussSurf(totalSurfQuadPoints, nodalStencil.size());
    Eigen::MatrixXd gaussSurfCoords(totalSurfQuadPoints, 2);
    Eigen::VectorXd gaussSurfWeights(totalSurfQuadPoints);

    // refCoord[0,1,2] are fixed to an arbitrary node to construct 2d integration matrix
    double refCoord[5];
    refCoord[0] = -1;
    refCoord[1] = -1;
    refCoord[2] = -1;
    std::vector<double> basisAtPoint(nodalStencil.size());
    double weightScale = 0.5*grid.getDx(3)*0.5*grid.getDx(4);
    for (int gaussIndexOuter = 0; gaussIndexOuter < numGaussPoints1d; gaussIndexOuter++)
    {
      refCoord[3] = gaussPoints1d[gaussIndexOuter];
      for (int gaussIndexInner = 0; gaussIndexInner < numGaussPoints1d; gaussIndexInner++)
      {
        refCoord[4] = gaussPoints1d[gaussIndexInner];
        int linIndex = gaussIndexOuter*numGaussPoints1d + gaussIndexInner;
        // Store coordinate of quadrature point (on ref. element)
        gaussSurfCoords.row(linIndex) << refCoord[3], refCoord[4];
        // Store integration weight of quadrature point (real space)
        gaussSurfWeights(linIndex) = weightScale*gaussWeights1d[gaussIndexOuter]*gaussWeights1d[gaussIndexInner];
        // Evaluate relevant basis functions at quadrature point
        nodalBasis5d->evalBasis(refCoord, basisAtPoint, nodalStencil);
        // Store results in matrix
        for (int nodeIndex = 0; nodeIndex < basisAtPoint.size(); nodeIndex++)
          gaussSurf(linIndex, nodeIndex) = basisAtPoint[nodeIndex];
      }
    }
 
    // Compute integration matrices (zeroth and first moments)
    // First element of output is the zeroth moment over the cell
    // Second element of output is the first v_para moment over the cell
    momentMatrix = Eigen::MatrixXd(2, nodalStencil.size());

    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
    {
      double mom0Val = 0.0;
      double mom1Val = 0.0;
      for (int quadIndex = 0; quadIndex < totalSurfQuadPoints; quadIndex++)
      {
        double vVal = 0.5*grid.getDx(3)*gaussSurfCoords(quadIndex,0);
        mom0Val += gaussSurfWeights(quadIndex)*gaussSurf(quadIndex, nodeIndex);
        mom1Val += gaussSurfWeights(quadIndex)*gaussSurf(quadIndex, nodeIndex)*vVal;
      }
      momentMatrix(0,nodeIndex) = mom0Val;
      momentMatrix(1,nodeIndex) = mom1Val;
    }

    // Fill out the node numbers on lower and upper surfaces in z
    lowerEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfLowerNodes(2));
    upperEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfUpperNodes(2));
    nodalBasis3d->getSurfLowerNodeNums(2, lowerEdgeNodeNums);
    nodalBasis3d->getSurfUpperNodeNums(2, upperEdgeNodeNums);
  }

  Lucee::UpdaterStatus
  MomentsAtEdges5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // Output parallel velocity moments (0, 1) vs time
    Lucee::Field<3, double>& outputMoment = this->getOut<Lucee::Field<3, double> >(0);
    
    outputMoment= 0.0; // clear out current contents

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr(); // for skin-cell
    Lucee::FieldPtr<double> outputMomentPtr = outputMoment.createPtr(); // Output pointer

    unsigned numNodes = nodalBasis5d->getNumNodes();

    double cellCentroid[5];
    int idx[5];

    // Check to see if we should integrate on the lower z plane
    if (localRgn.getLower(2) == globalRgn.getLower(2))
    {
      // Create a sequencer to loop over (x,y,z=0,v,mu) plane
      Lucee::RowMajorSequencer<5> seqLowerDim(localRgn.deflate(2));
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper value
        idx[2] = localRgn.getLower(2);
        distfIn.setPtr(distfInPtr, idx);
        outputMoment.setPtr(outputMomentPtr, idx[0], idx[1], idx[2]);
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Loop over the four configuration space vertices in this cell (specific to linear elements)
        Eigen::VectorXd distfReduced(nodalStencil.size());
        // V_PARA < 0.0 and on lower edge, so do the skin cell contribution
        if (cellCentroid[3] < 0.0)
        {
          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = lowerEdgeNodeNums[configNode];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            // Compute zeroth and first parallel velocity moments in this cell
            Eigen::VectorXd momentVector = momentMatrix*distfReduced;
            // Accumulate results (v = v_c*<1> + <v'>)
            outputMomentPtr[ lowerEdgeNodeNums[configNode] ] += cellCentroid[3]*momentVector(0) + momentVector(1);
          }
        }
        else if (integrateGhosts == true)
        {
          // V_PARA > 0.0 and on lower edge, so do the ghost cell contribution only if asked for
          idx[2] = localRgn.getLower(2)-1;
          distfIn.setPtr(distfInPtr, idx);
          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = upperEdgeNodeNums[configNode];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            // Compute zeroth and first parallel velocity moments in this cell
            Eigen::VectorXd momentVector = momentMatrix*distfReduced;
            // Accumulate results (v = v_c*<1> + <v'>)
            // Goes into lowerEdgeNodeNums since we are putting data only in skin cells
            outputMomentPtr[ lowerEdgeNodeNums[configNode] ] += cellCentroid[3]*momentVector(0) + momentVector(1);
          }
        }
      }
      // TESTING: Print out outputMomentPtr
      /*for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
        for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
        {
          outputMoment.setPtr(outputMomentPtr, ix, iy, globalRgn.getLower(2));
          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
            std::cout << "config " << lowerEdgeNodeNums[configNode] << " = " << outputMomentPtr[ lowerEdgeNodeNums[configNode] ] << std::endl;
        }*/
    }
    
    // Same as above loop, but done for upper edge in z
    if (localRgn.getUpper(2) == globalRgn.getUpper(2))
    {
      // Create a sequencer to loop over (x,y,z=Z_UPPER,v,mu) plane
      Lucee::RowMajorSequencer<5> seqLowerDim(localRgn.deflate(2));
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper valuea
        // the "-1" is because getUpper() gives index one beyond last domain cell
        idx[2] = localRgn.getUpper(2)-1;
        distfIn.setPtr(distfInPtr, idx);
        outputMoment.setPtr(outputMomentPtr, idx[0], idx[1], idx[2]);
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Loop over the four configuration space vertices in this cell (specific to linear elements)
        Eigen::VectorXd distfReduced(nodalStencil.size());
        // V_PARA > 0.0 and on upper edge, so do the skin cell contribution
        if (cellCentroid[3] > 0.0)
        {
          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = upperEdgeNodeNums[configNode];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            // Compute zeroth and first parallel velocity moments in this cell
            Eigen::VectorXd momentVector = momentMatrix*distfReduced;
            // Accumulate results (v = v_c*<1> + <v'>)
            outputMomentPtr[ upperEdgeNodeNums[configNode] ] += cellCentroid[3]*momentVector(0) + momentVector(1);
            double cellIncrement = cellCentroid[3]*momentVector(0) + momentVector(1);
            if (cellIncrement < 0.0)
            {
              for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              std::cout << "distfReduced = " << distfReduced << std::endl;
              std::cout << "cellIncrement = " << cellIncrement << std::endl;
              std::cout << "cell " << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << "," << idx[4] << std::endl;
              std::cout << "outputMoment[" << upperEdgeNodeNums[configNode] << "] new val = " << outputMomentPtr[ upperEdgeNodeNums[configNode] ] << std::endl;
            }
          }
        }
        else if (integrateGhosts == true)
        {
          // V_PARA < 0.0 and on upper edge, so do the ghost cell contribution only if asked for
          idx[2] = localRgn.getUpper(2);
          distfIn.setPtr(distfInPtr, idx);
          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = lowerEdgeNodeNums[configNode];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            // Compute zeroth and first parallel velocity moments in this cell
            Eigen::VectorXd momentVector = momentMatrix*distfReduced;
            // Accumulate results (v = v_c*<1> + <v'>)
            // Goes into upperEdgeNodeNums since we are putting data only in skin cells
            outputMomentPtr[ upperEdgeNodeNums[configNode] ] += cellCentroid[3]*momentVector(0) + momentVector(1);
          }
        }
      }
      // TESTING: Print out outputMomentPtr
      /*for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
        for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
        {
          outputMoment.setPtr(outputMomentPtr, ix, iy, globalRgn.getUpper(2)-1);
          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
            std::cout << "config " << upperEdgeNodeNums[configNode] << " = " << outputMomentPtr[ upperEdgeNodeNums[configNode] ] << std::endl;
        }
      */
    }
    return Lucee::UpdaterStatus();
  }

  void
  MomentsAtEdges5DUpdater::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Output: Degree 1 parallel velocity moment valid only on
    // upper and lower planes in z
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  MomentsAtEdges5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

  bool
  MomentsAtEdges5DUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

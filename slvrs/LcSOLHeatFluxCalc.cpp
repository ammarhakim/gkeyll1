/**
 * @file	LcSOLHeatFluxCalc.cpp
 *
 * @brief	Computes heat flux both edges as 2d fields for 5d sol simulations
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLHeatFluxCalc.h>
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
  const char *SOLHeatFluxCalc::id = "SOLHeatFluxCalc";

  SOLHeatFluxCalc::SOLHeatFluxCalc()
  {
  }

  void
  SOLHeatFluxCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLHeatFluxCalc::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLHeatFluxCalc::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLHeatFluxCalc::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLHeatFluxCalc::initialize()
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
 
    // Get (v,mu) quadrature data from 2d element
    int nSurfQuad = nodalBasis2d->getNumGaussNodes();
    int nlocal2d = nodalBasis2d->getNumNodes();
    interpMatrix2d = Eigen::MatrixXd(nSurfQuad, nlocal2d);
    gaussWeights2d = std::vector<double>(nSurfQuad);
    Lucee::Matrix<double> tempVolCoords(nSurfQuad, 3);
    Lucee::Matrix<double> tempVolQuad(nSurfQuad, nlocal2d);

    nodalBasis2d->getGaussQuadData(tempVolQuad, tempVolCoords, gaussWeights2d);
    copyLuceeToEigen(tempVolQuad, interpMatrix2d);

    // Scale gaussWeights2d to the right values since grid is (v,mu), not (x,y)
    double scaleCorrection = grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));
    for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
      gaussWeights2d[quadIndex] = scaleCorrection*gaussWeights2d[quadIndex];

    // Fill out the node numbers on lower and upper surfaces in z
    lowerEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfLowerNodes(2));
    upperEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfUpperNodes(2));
    nodalBasis3d->getSurfLowerNodeNums(2, lowerEdgeNodeNums);
    nodalBasis3d->getSurfUpperNodeNums(2, upperEdgeNodeNums);

    // Indicates if we will integrate in ghost cells in position space
    integrateGhosts = true;
  }

  Lucee::UpdaterStatus
  SOLHeatFluxCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // Hamiltonian parallel velocity derivative
    const Lucee::Field<5, double>& hamilDerivIn = this->getInp<Lucee::Field<5, double> >(1);
    // Hamiltonian
    const Lucee::Field<5, double>& hamilIn = this->getInp<Lucee::Field<5, double> >(2);
    // 3d magnetic field
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(3);
    // Heat flux on lower and upper edges as 2d fields
    Lucee::Field<2, double>& heatFluxLower = this->getOut<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& heatFluxUpper = this->getOut<Lucee::Field<2, double> >(1);
    
    heatFluxLower = 0.0; // clear out current contents
    heatFluxUpper = 0.0; // clear out current contents

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilDerivInPtr = hamilDerivIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilInPtr = hamilIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();

    Lucee::FieldPtr<double> heatFluxLowerPtr = heatFluxLower.createPtr(); // Output pointer
    Lucee::FieldPtr<double> heatFluxUpperPtr = heatFluxUpper.createPtr(); // Output pointer

    unsigned nlocal = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nSurfQuad = nodalBasis2d->getNumGaussNodes();

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
        heatFluxLower.setPtr(heatFluxLowerPtr, idx[0], idx[1]);
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Loop over the four configuration space vertices in this cell (specific to linear elements)
        Eigen::VectorXd distfReduced(nodalStencil.size());
        Eigen::VectorXd hamilDerivReduced(nodalStencil.size());
        Eigen::VectorXd hamilReduced(nodalStencil.size());
        Eigen::VectorXd distfReducedAtQuad(nSurfQuad);
        Eigen::VectorXd hamilDerivReducedAtQuad(nSurfQuad);
        Eigen::VectorXd hamilReducedAtQuad(nSurfQuad);
        // V_PARA < 0.0 and on lower edge, so do the skin cell contribution
        if (cellCentroid[3] < 0.0)
        {
          distfIn.setPtr(distfInPtr, idx);
          hamilDerivIn.setPtr(hamilDerivInPtr, idx);
          hamilIn.setPtr(hamilInPtr, idx);
          bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);

          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = lowerEdgeNodeNums[configNode];
            double bFieldVal = bFieldInPtr[configNodeIndex];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
            {
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilDerivReduced(nodeIndex) = hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilReduced(nodeIndex) = hamilInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            }

            // Compute fields at quadrature points
            distfReducedAtQuad = interpMatrix2d*distfReduced;
            hamilDerivReducedAtQuad = interpMatrix2d*hamilDerivReduced;
            hamilReducedAtQuad = interpMatrix2d*hamilReduced;
            // Compute 2d integral
            double integralResult = 0.0;
            for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
              integralResult += scaleFactor*bFieldVal*gaussWeights2d[quadIndex]*
                hamilReducedAtQuad(quadIndex)*hamilDerivReducedAtQuad(quadIndex)*distfReducedAtQuad(quadIndex);

            // Accumulate results of integration
            heatFluxLowerPtr[configNode] += integralResult;
          }
        }
        else if (integrateGhosts == true)
        {
          // V_PARA > 0.0 and on lower edge, so do the ghost cell contribution only if asked for
          idx[2] = localRgn.getLower(2)-1;
          distfIn.setPtr(distfInPtr, idx);
          hamilDerivIn.setPtr(hamilDerivInPtr, idx);
          hamilIn.setPtr(hamilInPtr, idx);
          bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);

          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = upperEdgeNodeNums[configNode];
            double bFieldVal = bFieldInPtr[configNodeIndex];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
            {
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilDerivReduced(nodeIndex) = hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilReduced(nodeIndex) = hamilInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            }

            // Compute fields at quadrature points
            distfReducedAtQuad = interpMatrix2d*distfReduced;
            hamilDerivReducedAtQuad = interpMatrix2d*hamilDerivReduced;
            hamilReducedAtQuad = interpMatrix2d*hamilReduced;
            // Compute 2d integral
            double integralResult = 0.0;
            for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
              integralResult += scaleFactor*bFieldVal*gaussWeights2d[quadIndex]*
                hamilReducedAtQuad(quadIndex)*hamilDerivReducedAtQuad(quadIndex)*distfReducedAtQuad(quadIndex);

            // Accumulate results of integration
            heatFluxLowerPtr[configNode] += integralResult;
          }
        }
      }
    }
    
    // Same as above loop, but done for upper edge in z
    if (localRgn.getUpper(2) == globalRgn.getUpper(2))
    {
      // Create a sequencer to loop over (x,y,z=Z_UPPER,v,mu) plane
      Lucee::RowMajorSequencer<5> seqLowerDim(localRgn.deflate(2));
      
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper value
        idx[2] = localRgn.getUpper(2)-1;
        heatFluxUpper.setPtr(heatFluxUpperPtr, idx[0], idx[1]);
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Loop over the four configuration space vertices in this cell (specific to linear elements)
        Eigen::VectorXd distfReduced(nodalStencil.size());
        Eigen::VectorXd hamilDerivReduced(nodalStencil.size());
        Eigen::VectorXd hamilReduced(nodalStencil.size());
        Eigen::VectorXd distfReducedAtQuad(nSurfQuad);
        Eigen::VectorXd hamilDerivReducedAtQuad(nSurfQuad);
        Eigen::VectorXd hamilReducedAtQuad(nSurfQuad);
        // V_PARA > 0.0 and on upper edge, so do the skin cell contribution
        if (cellCentroid[3] > 0.0)
        {
          distfIn.setPtr(distfInPtr, idx);
          hamilDerivIn.setPtr(hamilDerivInPtr, idx);
          hamilIn.setPtr(hamilInPtr, idx);
          bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);

          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = upperEdgeNodeNums[configNode];
            double bFieldVal = bFieldInPtr[configNodeIndex];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
            {
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilDerivReduced(nodeIndex) = hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilReduced(nodeIndex) = hamilInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            }

            // Compute fields at quadrature points
            distfReducedAtQuad = interpMatrix2d*distfReduced;
            hamilDerivReducedAtQuad = interpMatrix2d*hamilDerivReduced;
            hamilReducedAtQuad = interpMatrix2d*hamilReduced;
            // Compute 2d integral
            double integralResult = 0.0;
            for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
              integralResult += scaleFactor*bFieldVal*gaussWeights2d[quadIndex]*
                hamilReducedAtQuad(quadIndex)*hamilDerivReducedAtQuad(quadIndex)*distfReducedAtQuad(quadIndex);

            // Accumulate results of integration
            heatFluxUpperPtr[configNode] += integralResult;
          }
        }
        else if (integrateGhosts == true)
        {
          // V_PARA < 0.0 and on upper edge, so do the ghost cell contribution only if asked for
          idx[2] = localRgn.getUpper(2);
          distfIn.setPtr(distfInPtr, idx);
          hamilDerivIn.setPtr(hamilDerivInPtr, idx);
          hamilIn.setPtr(hamilInPtr, idx);
          bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);

          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = lowerEdgeNodeNums[configNode];
            double bFieldVal = bFieldInPtr[configNodeIndex];
            // At this particular configuration space vertix, copy all
            // nodes that occupy this location to a vector
            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
            {
              distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilDerivReduced(nodeIndex) = hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex];
              hamilReduced(nodeIndex) = hamilInPtr[nodalStencil[nodeIndex] + configNodeIndex];
            }

            // Compute fields at quadrature points
            distfReducedAtQuad = interpMatrix2d*distfReduced;
            hamilDerivReducedAtQuad = interpMatrix2d*hamilDerivReduced;
            hamilReducedAtQuad = interpMatrix2d*hamilReduced;
            // Compute 2d integral
            double integralResult = 0.0;
            for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
              integralResult += scaleFactor*bFieldVal*gaussWeights2d[quadIndex]*
                hamilReducedAtQuad(quadIndex)*hamilDerivReducedAtQuad(quadIndex)*distfReducedAtQuad(quadIndex);

            // Accumulate results of integration
            heatFluxUpperPtr[configNode] += integralResult;
          }
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLHeatFluxCalc::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: hamiltonian derivative
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: 3d magnetic field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: Heat flux on lower surface
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    // Output: Heat flux on upper surface
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  SOLHeatFluxCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLHeatFluxCalc::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

/**
 * @file	LcSOLUpperXPotentialBcUpdater.cpp
 *
 * @brief	Sets boundary conditions on potential for 5D SOL simulations
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLUpperXPotentialBcUpdater.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *SOLUpperXPotentialBcUpdater::id = "SOLUpperXPotentialBcUpdater";

  SOLUpperXPotentialBcUpdater::SOLUpperXPotentialBcUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  SOLUpperXPotentialBcUpdater::~SOLUpperXPotentialBcUpdater()
  {
  }
  
  void
  SOLUpperXPotentialBcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("SOLUpperXPotentialBcUpdater::readInput: Must specify 2D element to use using 'basis'");

    // by default, only set potential to a constant on the upper x edge
    // here is an option to also set potential to a constant lower x edge
    if (tbl.hasBool("applyLowerEdge"))
      applyLowerEdge = tbl.getBool("applyLowerEdge");
    else applyLowerEdge = false;
  }

  void
  SOLUpperXPotentialBcUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 2D
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    // get volume interpolation matrices for 2d element
    int nSurfQuad2d = nodalBasis2d->getNumSurfGaussNodes();
    
    std::vector<double> surfWeightsUpper2d(nSurfQuad2d);
    std::vector<double> surfWeightsLower2d(nSurfQuad2d);
    Lucee::Matrix<double> tempSurfQuad2d(nSurfQuad2d, nlocal2d);
    Lucee::Matrix<double> tempSurfCoords2d(nSurfQuad2d, 3);

    nodalBasis2d->getSurfUpperGaussQuadData(0, tempSurfQuad2d, tempSurfCoords2d, surfWeightsUpper2d);
    Eigen::MatrixXd surfInterpUpper2d(nSurfQuad2d, nlocal2d);
    copyLuceeToEigen(tempSurfQuad2d, surfInterpUpper2d);

    nodalBasis2d->getSurfLowerGaussQuadData(0, tempSurfQuad2d, tempSurfCoords2d, surfWeightsLower2d);
    Eigen::MatrixXd surfInterpLower2d(nSurfQuad2d, nlocal2d);
    copyLuceeToEigen(tempSurfQuad2d, surfInterpLower2d);
    
    // when multiplied by the solution in an element, gives
    // potential integrated on lower and upper sufaces in x
    integrationMatrixUpper = Eigen::VectorXd(nlocal2d);
    integrationMatrixLower = Eigen::VectorXd(nlocal2d);

    for (int i = 0; i < nlocal2d; i++)
    {
      // Each row of interpolation matrix is a different interpolation point
      // Each column is a different basis function
      double basisIntegralUpper = 0.0;
      double basisIntegralLower = 0.0;
      // Calculate integral on lower and upper surfaces
      for (int quadIndex = 0; quadIndex < surfInterpUpper2d.rows(); quadIndex++)
      {
        basisIntegralUpper += surfWeightsUpper2d[quadIndex]*surfInterpUpper2d(quadIndex, i);
        basisIntegralLower += surfWeightsLower2d[quadIndex]*surfInterpLower2d(quadIndex, i);
      }

      // Store result
      integrationMatrixUpper(i) = basisIntegralUpper;
      integrationMatrixLower(i) = basisIntegralLower;
    }

    upperEdgeNodeNums = std::vector<int>(nodalBasis2d->getNumSurfUpperNodes(0));
    nodalBasis2d->getSurfUpperNodeNums(0, upperEdgeNodeNums);

    lowerEdgeNodeNums = std::vector<int>(nodalBasis2d->getNumSurfLowerNodes(0));
    nodalBasis2d->getSurfLowerNodeNums(0, lowerEdgeNodeNums);
  }

  Lucee::UpdaterStatus
  SOLUpperXPotentialBcUpdater::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // get input field (2d)
    const Lucee::Field<2, double>& phiAvgIn = this->getInp<Lucee::Field<2, double> >(0);
    // get output field (2d)
    Lucee::Field<2, double>& phiLowerOut = this->getOut<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& phiUpperOut = this->getOut<Lucee::Field<2, double> >(1);
    // output dynvector
    Lucee::DynVector<double>& avgPhiAtEdge = this->getOut<Lucee::DynVector<double> >(2);

    // Local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
    // Global region to tell if we are on domain boundary
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    // iterators into fields
    Lucee::ConstFieldPtr<double> phiAvgPtr = phiAvgIn.createConstPtr();
    Lucee::FieldPtr<double> phiLowerPtr = phiLowerOut.createPtr();
    Lucee::FieldPtr<double> phiUpperPtr = phiUpperOut.createPtr();
 
    int idx[2];
    Lucee::RowMajorSequencer<2> seq(localRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    double localIntUpper = 0.0;

    // Calculate the average value of potential on upper x surface
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);

      if ( idx[0] == globalRgn.getUpper(0)-1 )
      {
        phiAvgIn.setPtr(phiAvgPtr, idx);

        Eigen::VectorXd phiAvgVec(nlocal2d);
        for (int i = 0; i < nlocal2d; i++)
          phiAvgVec(i) = phiAvgPtr[i];

        localIntUpper += phiAvgVec.dot(integrationMatrixUpper);
      }
    }

    double totalIntUpper = localIntUpper;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    // Every processor has a copy of the average potential
    comm->allreduce(1, &localIntUpper, &totalIntUpper, TX_SUM);

    // 4-16-16: To match torpex, set potential to zero
    totalIntUpper = totalIntUpper / (grid.getDx(1)*(globalRgn.getUpper(1)-globalRgn.getLower(1)));
    totalIntUpper = 0.0;

    seq.reset();

    // Only replace potential with average potential if it is an upper surface node
    while(seq.step())
    {
      seq.fillWithIndex(idx);

      if ( idx[0] == globalRgn.getUpper(0)-1 )
      {
        phiLowerOut.setPtr(phiLowerPtr, idx);
        phiUpperOut.setPtr(phiUpperPtr, idx);

        for (int nodeIndex = 0; nodeIndex < upperEdgeNodeNums.size(); nodeIndex++)
        {
          phiLowerPtr[upperEdgeNodeNums[nodeIndex]] = totalIntUpper;
          phiUpperPtr[upperEdgeNodeNums[nodeIndex]] = totalIntUpper;
        }
      }
    }

    // Do the same thing on the lower edge if needed
    if (applyLowerEdge == true)
    {
      double localIntLower = 0.0;
      // Calculate the average value of potential on lower x surface
      seq.reset();
      while(seq.step())
      {
        seq.fillWithIndex(idx);
        grid.setIndex(idx);

        if ( idx[0] == globalRgn.getLower(0) )
        {
          phiAvgIn.setPtr(phiAvgPtr, idx);

          Eigen::VectorXd phiAvgVec(nlocal2d);
          for (int i = 0; i < nlocal2d; i++)
            phiAvgVec(i) = phiAvgPtr[i];

          localIntLower += phiAvgVec.dot(integrationMatrixLower);
        }
      }

      double totalIntLower = localIntLower;
      // Every processor has a copy of the average potential
      comm->allreduce(1, &localIntLower, &totalIntLower, TX_SUM);

      // 4-16-16 to match TORPEX, set potential at edges to zero
      totalIntLower = 0.0;//totalIntLower / (grid.getDx(1)*(globalRgn.getUpper(1)-globalRgn.getLower(1)));

      seq.reset();

      // Only replace potential with average potential if it is an upper surface node
      while(seq.step())
      {
        seq.fillWithIndex(idx);

        if ( idx[0] == globalRgn.getLower(0) )
        {
          phiLowerOut.setPtr(phiLowerPtr, idx);
          phiUpperOut.setPtr(phiUpperPtr, idx);

          for (int nodeIndex = 0; nodeIndex < lowerEdgeNodeNums.size(); nodeIndex++)
          {
            phiLowerPtr[lowerEdgeNodeNums[nodeIndex]] = totalIntLower;
            phiUpperPtr[lowerEdgeNodeNums[nodeIndex]] = totalIntLower;
          }
        }
      }

      std::vector<double> writeVal(2);
      writeVal[0] = totalIntUpper;
      writeVal[1] = totalIntLower;
      avgPhiAtEdge.appendData(t, writeVal);
    }
    else
    {
      // just write a single value
      std::vector<double> writeVal(1);
      writeVal[0] = totalIntUpper;
      avgPhiAtEdge.appendData(t, writeVal);
    }

    // Average boundary nodes that are not on corners
    seq.reset();

    // Additional pointers needed for averaging operation
    Lucee::FieldPtr<double> phiLowerPtr_r = phiLowerOut.createPtr();
    Lucee::FieldPtr<double> phiUpperPtr_r = phiUpperOut.createPtr();
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      // Average lower 0 surface nodes only if applyLowerEdge == false
      if ((applyLowerEdge == false) && (idx[0] == globalRgn.getLower(0) && idx[1] < globalRgn.getUpper(1)-1 ))
      {
        if(idx[1] == localRgn.getLower(1) && localRgn.getLower(1) != globalRgn.getLower(1))
        {
          // Bottom left node of this region is not a corner
          // Should average this node with the ghost cell node below it
          phiLowerOut.setPtr(phiLowerPtr, idx[0], idx[1]-1);
          phiLowerOut.setPtr(phiLowerPtr_r, idx);
          double avgVal = 0.5*(phiLowerPtr[2] + phiLowerPtr_r[0]);
          phiLowerPtr[2] = avgVal;
          phiLowerPtr_r[0] = avgVal;
          // Duplicate code
          phiUpperOut.setPtr(phiUpperPtr, idx[0], idx[1]-1);
          phiUpperOut.setPtr(phiUpperPtr_r, idx);
          avgVal = 0.5*(phiUpperPtr[2] + phiUpperPtr_r[0]);
          phiUpperPtr[2] = avgVal;
          phiUpperPtr_r[0] = avgVal;
        }
        // Average upper left node of this cell with bottom left node of the cell above it
        phiLowerOut.setPtr(phiLowerPtr, idx);
        phiLowerOut.setPtr(phiLowerPtr_r, idx[0], idx[1]+1);
        double avgVal = 0.5*(phiLowerPtr[2] + phiLowerPtr_r[0]);
        phiLowerPtr[2] = avgVal;
        phiLowerPtr_r[0] = avgVal;
        // Duplicate code
        phiUpperOut.setPtr(phiUpperPtr, idx);
        phiUpperOut.setPtr(phiUpperPtr_r, idx[0], idx[1]+1);
        avgVal = 0.5*(phiUpperPtr[2] + phiUpperPtr_r[0]);
        phiUpperPtr[2] = avgVal;
        phiUpperPtr_r[0] = avgVal;
      }
      // Average lower 1 surface nodes
      if (idx[1] == globalRgn.getLower(1) && idx[0] < globalRgn.getUpper(0)-1)
      {
        if(idx[0] == localRgn.getLower(0) && localRgn.getLower(0) != globalRgn.getLower(0))
        {
          // Bottom left node of this region is not a corner
          // Should average this node with the ghost cell node to the left of it
          phiLowerOut.setPtr(phiLowerPtr, idx[0]-1, idx[1]);
          phiLowerOut.setPtr(phiLowerPtr_r, idx);
          double avgVal = 0.5*(phiLowerPtr[1] + phiLowerPtr_r[0]);
          phiLowerPtr[1] = avgVal;
          phiLowerPtr_r[0] = avgVal;
          // Duplicate code
          phiUpperOut.setPtr(phiUpperPtr, idx[0]-1, idx[1]);
          phiUpperOut.setPtr(phiUpperPtr_r, idx);
          avgVal = 0.5*(phiUpperPtr[1] + phiUpperPtr_r[0]);
          phiUpperPtr[1] = avgVal;
          phiUpperPtr_r[0] = avgVal;
        }
        // Average bottom right node with bottom left node of adjacent cell
        phiLowerOut.setPtr(phiLowerPtr, idx);
        phiLowerOut.setPtr(phiLowerPtr_r, idx[0]+1, idx[1]);
        double avgVal = 0.5*(phiLowerPtr[1] + phiLowerPtr_r[0]);
        phiLowerPtr[1] = avgVal;
        phiLowerPtr_r[0] = avgVal;
        // Duplicate code
        phiUpperOut.setPtr(phiUpperPtr, idx);
        phiUpperOut.setPtr(phiUpperPtr_r, idx[0]+1, idx[1]);
        avgVal = 0.5*(phiUpperPtr[1] + phiUpperPtr_r[0]);
        phiUpperPtr[1] = avgVal;
        phiUpperPtr_r[0] = avgVal;
      }
      // Average upper 1 surface nodes
      if (idx[1] == globalRgn.getUpper(1)-1 && idx[0] < globalRgn.getUpper(0)-1)
      {
        if(idx[0] == localRgn.getLower(0) && localRgn.getLower(0) != globalRgn.getLower(0))
        {
          // Bottom top left node of this region is not a corner
          // Should average this node with the ghost cell node to the left of it
          phiLowerOut.setPtr(phiLowerPtr, idx[0]-1, idx[1]);
          phiLowerOut.setPtr(phiLowerPtr_r, idx);
          double avgVal = 0.5*(phiLowerPtr[3] + phiLowerPtr_r[2]);
          phiLowerPtr[3] = avgVal;
          phiLowerPtr_r[2] = avgVal;
          // Duplicate code
          phiUpperOut.setPtr(phiUpperPtr, idx[0]-1, idx[1]);
          phiUpperOut.setPtr(phiUpperPtr_r, idx);
          avgVal = 0.5*(phiUpperPtr[3] + phiUpperPtr_r[2]);
          phiUpperPtr[3] = avgVal;
          phiUpperPtr_r[2] = avgVal;
        }
        // Average upper right node with upper left node of adjacent cell
        phiLowerOut.setPtr(phiLowerPtr, idx);
        phiLowerOut.setPtr(phiLowerPtr_r, idx[0]+1, idx[1]);
        double avgVal = 0.5*(phiLowerPtr[3] + phiLowerPtr_r[2]);
        phiLowerPtr[3] = avgVal;
        phiLowerPtr_r[2] = avgVal;
        // Duplicate code
        phiUpperOut.setPtr(phiUpperPtr, idx);
        phiUpperOut.setPtr(phiUpperPtr_r, idx[0]+1, idx[1]);
        avgVal = 0.5*(phiUpperPtr[3] + phiUpperPtr_r[2]);
        phiUpperPtr[3] = avgVal;
        phiUpperPtr_r[2] = avgVal;
      }
    }
    
    return Lucee::UpdaterStatus();
  }

  void
  SOLUpperXPotentialBcUpdater::declareTypes()
  {
    // Input potential on a 2d (x,y) field
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Zonal-averaged potential on a 2d (x,z) field
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    // Value of zonal-averaged sheath potential
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  SOLUpperXPotentialBcUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}

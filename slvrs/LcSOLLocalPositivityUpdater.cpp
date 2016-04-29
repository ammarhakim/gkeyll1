/**
 * @file	LcSOLLocalPositivityUpdater.cpp
 *
 * @brief	Updater to enforce positivity preservation for 5d SOL.
 * Only reallocates solution across (v,mu) components in a cell.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLLocalPositivityUpdater.h>
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
  const char *SOLLocalPositivityUpdater::id = "SOLLocalPositivityUpdater";

  SOLLocalPositivityUpdater::SOLLocalPositivityUpdater()
  {
  }

  void
  SOLLocalPositivityUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLLocalPositivityUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLLocalPositivityUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLLocalPositivityUpdater::readInput: Must specify element to use using 'basis2d'");
  }

  void
  SOLLocalPositivityUpdater::initialize()
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

    // Construct vector to integrate (v,mu) cell
    mom0Vector = Eigen::VectorXd(nlocal2d);
    // Integrate each basis function over entire cell
    for (int nodeIndex = 0; nodeIndex < nlocal2d; nodeIndex++)
    {
      double integralResult = 0.0;
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        integralResult += gaussWeights2d[quadIndex]*interpMatrix2d(quadIndex, nodeIndex);
      mom0Vector(nodeIndex) = integralResult;
    }
  }

  Lucee::UpdaterStatus
  SOLLocalPositivityUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function to be modified
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr();

    unsigned nlocal = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nSurfQuad = nodalBasis2d->getNumGaussNodes();

    double cellCentroid[5];
    int idx[5];

    // Loop over all local cells
    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      // Loop over the four configuration space vertices in this cell (specific to linear elements)
      Eigen::VectorXd distfReduced(nodalStencil.size());
      distfOut.setPtr(distfOutPtr, idx);

      for (int configNode = 0; configNode < nlocal3d; configNode++)
      {
        // At this particular configuration space vertix, copy all
        // nodes that occupy this location to a vector
        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
          distfReduced(nodeIndex) = distfOutPtr[nodalStencil[nodeIndex] + configNode];

        // Compute 2d integral for density before applying positivity
        double originalNum = mom0Vector.dot(distfReduced);

        if (originalNum < 0.0)
        {
          std::cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << "," << idx[4] 
            << ") originalNum is less than zero! " << originalNum << std::endl;
          std::cout << distfReduced << std::endl;
          return Lucee::UpdaterStatus(false, 0.0);
        }

        // Zero out distfVector entries that are negative
        for (int i = 0; i < distfReduced.size(); i++)
        {
          if (distfReduced(i) < 0.0)
            distfReduced(i) = 0.0;
        }

        // Compute 2d integral for density after applying positivity
        double modifiedNum = mom0Vector.dot(distfReduced);

        // Need to develop a check to see under what conditions it's okay to set modifiedNum to zero
        if (modifiedNum < 0.0)
        {
          // Set modifiedNum to zero if it won't change the density too much
          if (std::fabs(modifiedNum-originalNum) < 1e-8*std::fabs(modifiedNum))
            modifiedNum = 0.0;
          else
          {
            std::cout << "modifiedNum is negative and unable to set modifiedNum to 0." << std::endl;
            return Lucee::UpdaterStatus(false, 0.0);
          }
        }

        // Assume if originalNum is negative, it's very small so make zero
        if (originalNum < 0.0)
        {
          // Write modified values to distfOut
          for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
            distfOutPtr[nodalStencil[nodeIndex] + configNode] = std::fabs((originalNum/modifiedNum))*distfReduced(nodeIndex);
          //return Lucee::UpdaterStatus(false, 0.0);
        }
        else
        {
          // Write modified values to distfOut
          for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
            distfOutPtr[nodalStencil[nodeIndex] + configNode] = (originalNum/modifiedNum)*distfReduced(nodeIndex);
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLLocalPositivityUpdater::declareTypes()
  {
    // Output: distribution function to be modified
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLLocalPositivityUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLLocalPositivityUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

/**
 * @file	LcSOLPositivityScaleNodeUpdater.cpp
 *
 * @brief	Updater used to adjust energy at each node by use of a drag term
 * Uses result of LcSOLPositivityDragNodeUpdater to add the correct portion to match
 * the desired energy at each node.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLPositivityScaleNodeUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLPositivityScaleNodeUpdater::id = "SOLPositivityScaleNodeUpdater";

  SOLPositivityScaleNodeUpdater::SOLPositivityScaleNodeUpdater()
  {
  }

  void
  SOLPositivityScaleNodeUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLPositivityScaleNodeUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLPositivityScaleNodeUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLPositivityScaleNodeUpdater::readInput: Must specify element to use using 'basis2d'");
  }

  void
  SOLPositivityScaleNodeUpdater::initialize()
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
  SOLPositivityScaleNodeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function difference between positivity and drag+positivity steps
    const Lucee::Field<5, double>& distfDelta = this->getInp<Lucee::Field<5, double> >(0);
    // Energy at nodes before positivity (target energy)
    const Lucee::Field<3, double>& energyOrigIn = this->getInp<Lucee::Field<3, double> >(1);
    // Energy at nodes after positivity
    const Lucee::Field<3, double>& energyPosIn = this->getInp<Lucee::Field<3, double> >(2);
    // Energy at nodes after maximum drag
    const Lucee::Field<3, double>& energyDragIn = this->getInp<Lucee::Field<3, double> >(3);
    // Distribution function with positivity and correct drag applied
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfDeltaPtr   = distfDelta.createConstPtr();
    Lucee::ConstFieldPtr<double> energyOrigInPtr = energyOrigIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyPosInPtr  = energyPosIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyDragInPtr = energyDragIn.createConstPtr();

    // Remember not to clear distfOut
    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr(); // Output pointer

    unsigned nlocal = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nSurfQuad = nodalBasis2d->getNumGaussNodes();

    double cellCentroid[5];
    int idx[5];

    // Loop over each cell in configuration space
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      idx[0] = ix;
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        idx[1] = iy;
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[2] = iz;
          energyOrigIn.setPtr(energyOrigInPtr, idx[0], idx[1], idx[2]);
          energyPosIn.setPtr(energyPosInPtr, idx[0], idx[1], idx[2]);
          energyDragIn.setPtr(energyDragInPtr, idx[0], idx[1], idx[2]);

          // Loop over each node in configuration space
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // At this particular configuration node, loop in vPara and mu space
            // to scale distribution function only if energy at point needs to be corrected
            if (std::fabs(energyPosInPtr[configNode] - energyOrigInPtr[configNode]) > 
                1e-10*energyOrigInPtr[configNode])
            {
              // Compute scale factor to add delta f at this point
              double scaleFactor = (energyOrigInPtr[configNode]-energyPosInPtr[configNode])/
                (energyDragInPtr[configNode]-energyPosInPtr[configNode]);

              // Make sure scaleFactor is less than one
              if (scaleFactor > 1.0)
              {
                std::cout << "Drag term is larger than possible = " << scaleFactor << std::endl;
                scaleFactor = 1.0;
              }

              for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
              {
                idx[3] = iv;
                for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
                {
                  idx[4] = imu;

                  distfDelta.setPtr(distfDeltaPtr, idx);
                  distfOut.setPtr(distfOutPtr, idx);
                  // Set fOut to give the right energy at this node
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  {
                    distfOutPtr[nodalStencil[nodeIndex] + configNode] =
                      distfOutPtr[nodalStencil[nodeIndex] + configNode] +
                      scaleFactor*distfDeltaPtr[nodalStencil[nodeIndex] + configNode];
                    if (distfOutPtr[nodalStencil[nodeIndex] + configNode] < 0.0)
                      std::cout << "distfOutPtr = " << distfOutPtr[nodalStencil[nodeIndex] + configNode] << std::endl;
                  }
                }
              }
            }
          }
        }
      }
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLPositivityScaleNodeUpdater::declareTypes()
  {
    // Input: difference in distribution functions of positivity and drag steps
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: target energy of distribution function at nodes
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of distribition function at nodes after positivity
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of distribition function at nodes after maximum drag
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: non-negative distribution function with correct energy
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLPositivityScaleNodeUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLPositivityScaleNodeUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}
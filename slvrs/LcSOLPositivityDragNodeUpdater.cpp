/**
 * @file	LcSOLPositivityDragNodeUpdater.cpp
 *
 * @brief	Updater used to adjust energy at each node by use of a drag term
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLPositivityDragNodeUpdater.h>
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
  const char *SOLPositivityDragNodeUpdater::id = "SOLPositivityDragNodeUpdater";

  SOLPositivityDragNodeUpdater::SOLPositivityDragNodeUpdater()
  {
  }

  void
  SOLPositivityDragNodeUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLPositivityDragNodeUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLPositivityDragNodeUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLPositivityDragNodeUpdater::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLPositivityDragNodeUpdater::initialize()
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
  SOLPositivityDragNodeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function with incorrect energy
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // Energy of distfIn at each node
    const Lucee::Field<3, double>& energyModIn = this->getInp<Lucee::Field<3, double> >(1);
    // Desired energy at each node
    const Lucee::Field<3, double>& energyOrigIn = this->getInp<Lucee::Field<3, double> >(2);
    // Distribution function after drag term
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr       = distfIn.createConstPtr();
    // Another pointer to identify upwind neighbor
    Lucee::ConstFieldPtr<double> distfInUpwindPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyModInPtr   = energyModIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyOrigInPtr  = energyOrigIn.createConstPtr();

    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr(); // Output pointer
    distfOut = 0.0; // clear out current contents

    unsigned nlocal = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nSurfQuad = nodalBasis2d->getNumGaussNodes();

    double cellCentroid[5];
    int idx[5];
    int idxUpwind[5];
    Eigen::VectorXd distfReduced(nodalStencil.size());
    Eigen::VectorXd distfUpwindReduced(nodalStencil.size());

    double dt = t-this->getCurrTime();
    // Figure out value of parallel velocity on last edge that flux is computed on
    idx[0] = globalRgn.getUpper(0)-1;
    idx[1] = globalRgn.getUpper(1)-1;
    idx[2] = globalRgn.getUpper(2)-1;
    idx[3] = globalRgn.getUpper(3)-1;
    idx[4] = globalRgn.getUpper(4)-1;
    // Set grid to last parallel velocity cell (upper boundary)
    grid.setIndex(idx);
    grid.getCentroid(cellCentroid);
    double alpha = 1/(cellCentroid[3] - 0.5*grid.getDx(3));

    // Loop over each cell in configuration space
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      idx[0] = ix;
      idxUpwind[0] = ix;
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        idx[1] = iy;
        idxUpwind[1] = iy;
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[2] = iz;
          idxUpwind[2] = iz;
          energyModIn.setPtr(energyModInPtr, idx[0], idx[1], idx[2]);
          energyOrigIn.setPtr(energyOrigInPtr, idx[0], idx[1], idx[2]);
          // Loop over each node in configuration space
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // At this particular configuration node, loop in vPara and mu space
            // to perform drag step only if energy at point needs to be corrected
            if (std::fabs(energyModInPtr[configNode] - energyOrigInPtr[configNode]) > 
                1e-10*energyOrigInPtr[configNode])
            {
              // Temporary Diagnostic: calculate total density before and after drag:
              /*double numDensityBefore = 0.0;
              for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
              {
                idx[3] = iv;
                for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
                {
                  idx[4] = imu;
                  distfIn.setPtr(distfInPtr, idx);
                  // Construct reduced distf
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNode];
                  // Compute density of (v,mu) cell
                  numDensityBefore += mom0Vector.dot(distfReduced);
                }
              }*/

              for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
              {
                idx[3] = iv;
                for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
                {
                  idx[4] = imu;

                  // Get the coordinates of cell center
                  grid.setIndex(idx);
                  grid.getCentroid(cellCentroid);

                  distfIn.setPtr(distfInPtr, idx);
                  // Construct reduced distf
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNode];
                  // Compute density of (v,mu) cell
                  double fOldAvg = mom0Vector.dot(distfReduced);
                  double fLowerAvg = 0.0;
                  double fUpperAvg = 0.0;

                  if (fOldAvg == 0.0)
                  {
                    //std::cout << "fOldAvg is zero. This should probably not happen." << std::endl;
                    //std::cout << "distfReduced" << std::endl << distfReduced << std::endl;
                    continue;
                  }
                  else if (fOldAvg < 0.0)
                  {
                    std::cout << "fOldAvg is less than zero. This should probably not happen." << std::endl;
                    //std::cout << "distfReduced" << std::endl << distfReduced << std::endl;
                  }
                  
                  // Figure out which neighboring cell to use to calculate upwindAvg
                  if (cellCentroid[3] > 0.0)
                  {
                    fLowerAvg = fOldAvg;
                    if (iv < globalRgn.getUpper(3)-1)
                    {
                      idxUpwind[3] = idx[3] + 1;
                      idxUpwind[4] = idx[4];
                      distfIn.setPtr(distfInUpwindPtr, idxUpwind);
                      for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                        distfUpwindReduced(nodeIndex) = distfInUpwindPtr[nodalStencil[nodeIndex] + configNode];
                      fUpperAvg = mom0Vector.dot(distfUpwindReduced);
                    }
                  }
                  else if (cellCentroid[3] < 0.0)
                  {
                    fUpperAvg = fOldAvg;
                    if (iv > globalRgn.getLower(3))
                    {
                      idxUpwind[3] = idx[3] - 1;
                      idxUpwind[4] = idx[4];
                      distfIn.setPtr(distfInUpwindPtr, idxUpwind);
                      for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                        distfUpwindReduced(nodeIndex) = distfInUpwindPtr[nodalStencil[nodeIndex] + configNode];
                      fLowerAvg = mom0Vector.dot(distfUpwindReduced);
                    }
                  }

                  // Put in limiters on cell average?
                  if (fUpperAvg < 0.0)
                    fUpperAvg = 0.0;
                  if (fLowerAvg < 0.0)
                    fLowerAvg = 0.0;

                  double fIncrement = (cellCentroid[3] + 0.5*grid.getDx(3))*fUpperAvg -
                    (cellCentroid[3] - 0.5*grid.getDx(3))*fLowerAvg;
                  double fNewAvg = fOldAvg + alpha*fIncrement;

                  if (fNewAvg < 0.0)
                  {
                    // Set fNewAvg to 0 if negative (presumably it happens when answer is approximately zero)
                    if (std::fabs(fNewAvg/fOldAvg) < 1e-8)
                    {
                      fNewAvg = 0.0;
                    }
                    else
                    {
                    std::cout << "fNewAvg is negative! = " << fNewAvg << std::endl;
                    std::cout << "fOldAvg = " << fOldAvg << std::endl;
                    std::cout << "alpha*fIncrement = " << alpha*fIncrement << std::endl;
                    std::cout << "fIncrement upper = " << alpha*(cellCentroid[3] + 0.5*grid.getDx(3))*fUpperAvg << std::endl;
                    std::cout << "fIncrement lower = " << alpha*(cellCentroid[3] - 0.5*grid.getDx(3))*fLowerAvg << std::endl;
                    }
                  }

                  
                  distfOut.setPtr(distfOutPtr, idx);
                  // Set fOut by scaling fIn to have the right density
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    distfOutPtr[nodalStencil[nodeIndex] + configNode] = 
                      (fNewAvg/fOldAvg)*distfInPtr[nodalStencil[nodeIndex] + configNode];
                }
              }

              // Temporary Diagnostic: calculate total density before and after drag:
              /*double numDensityAfter = 0.0;
              for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
              {
                idx[3] = iv;
                for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
                {
                  idx[4] = imu;
                  distfOut.setPtr(distfOutPtr, idx);
                  // Construct reduced distf
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    distfReduced(nodeIndex) = distfOutPtr[nodalStencil[nodeIndex] + configNode];
                  // Compute density of (v,mu) cell
                  numDensityAfter += mom0Vector.dot(distfReduced);
                }
              }*/

              //std::cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << ") node " << configNode << std::endl;
              //std::cout << "numDensityBefore = " << numDensityBefore << std::endl;
              //std::cout << "numDensityAfter = " << numDensityAfter << std::endl;

            }
          }
        }
      }
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLPositivityDragNodeUpdater::declareTypes()
  {
    // Input: modified distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: energy of modified distribution function at nodes
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of original distribition function at nodes
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: another modified distribution function
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLPositivityDragNodeUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLPositivityDragNodeUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

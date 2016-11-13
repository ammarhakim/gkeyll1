/**
 * @file	LcSOLPositivityDragSliceInMuUpdater.cpp
 *
 * @brief	Updater used to adjust energy at each node by use of a drag term that acts on a 1d slice
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLPositivityDragSliceInMuUpdater.h>
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
  const char *SOLPositivityDragSliceInMuUpdater::id = "SOLPositivityDragSliceInMuUpdater";

  SOLPositivityDragSliceInMuUpdater::SOLPositivityDragSliceInMuUpdater()
  {
  }

  void
  SOLPositivityDragSliceInMuUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLPositivityDragSliceInMuUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLPositivityDragSliceInMuUpdater::readInput: Must specify element to use using 'basis3d'");
  }

  void
  SOLPositivityDragSliceInMuUpdater::initialize()
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

    // Should be 0,16
    std::vector<int> fixedVParDimList;
    fixedVParDimList.push_back(0);
    fixedVParDimList.push_back(1);
    fixedVParDimList.push_back(2);
    fixedVParDimList.push_back(3);
    dxMin = fixedVParDimList[0];
    for (int d = 1; d < fixedVParDimList.size(); d++)
      dxMin = std::min(dxMin, grid.getDx(fixedVParDimList[d]));
    // Find all nodes that share the same set of coordinates as node 0
    nodalStencilFixedVPar = std::vector<int>(nlocal);
    stencilIndex = 0;
    for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
    {
      if (sameSpecifiedCoords(0, nodeIndex, dxMin, nodeCoords, fixedVParDimList) == true)
      {
        nodalStencilFixedVPar[stencilIndex] = nodeIndex;
        stencilIndex++;
      }
    }
    nodalStencilFixedVPar.resize(stencilIndex);

    // Should be 0,8
    std::vector<int> fixedMuDimList;
    fixedMuDimList.push_back(0);
    fixedMuDimList.push_back(1);
    fixedMuDimList.push_back(2);
    fixedMuDimList.push_back(4);
    dxMin = fixedMuDimList[0];
    for (int d = 1; d < fixedMuDimList.size(); d++)
      dxMin = std::min(dxMin, grid.getDx(fixedMuDimList[d]));
    // Find all nodes that share the same set of coordinates as node 0
    nodalStencilFixedMu = std::vector<int>(nlocal);
    stencilIndex = 0;
    for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
    {
      if (sameSpecifiedCoords(0, nodeIndex, dxMin, nodeCoords, fixedMuDimList) == true)
      {
        nodalStencilFixedMu[stencilIndex] = nodeIndex;
        stencilIndex++;
      }
    }
    nodalStencilFixedMu.resize(stencilIndex);
  }

  Lucee::UpdaterStatus
  SOLPositivityDragSliceInMuUpdater::update(double t)
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

    double cellCentroid[5];
    int idx[5];
    int idxUpwind[5];
    Eigen::VectorXd distfReduced(nodalStencilFixedVPar.size());
    Eigen::VectorXd distfUpwindReduced(nodalStencilFixedVPar.size());

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
    // Zero flux boundary conditions, so use max velocity - dmu
    double maxPosMu = cellCentroid[4] - 0.5*grid.getDx(4);

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
                1e-12*energyOrigInPtr[configNode])
            {

              for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
              {
                idx[3] = iv;
                idxUpwind[3] = iv;
                for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
                {
                  idx[4] = imu;

                  // Get the coordinates of cell center
                  grid.setIndex(idx);
                  grid.getCentroid(cellCentroid);

                  distfIn.setPtr(distfInPtr, idx);

                  // Loop over each slice at fixed vPar
                  for (int vSliceIndex = 0; vSliceIndex < nodalStencilFixedMu.size(); vSliceIndex++)
                  {
                    // Construct reduced distf
                    for (int nodeIndex = 0; nodeIndex < nodalStencilFixedVPar.size(); nodeIndex++)
                      distfReduced(nodeIndex) = distfInPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[nodeIndex]];

                    // Compute density of distfReduced (for now just use average)
                    double fOldAvg = distfReduced.mean();

                    // Add density to lowest-energy cell only
                    idx[4] = localRgn.getLower(4);
                    distfOut.setPtr(distfOutPtr, idx);
                    distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[0]] = distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[0]] + 2*fOldAvg;
                    distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[1]] = 0.0;
                    
                    /*
                    double fLowerAvg = 0.0;
                    double fUpperAvg = 0.0;

                    if (fOldAvg < 0.0)
                    {
                      std::cout << "fOldAvg is less than zero. This should probably not happen." << std::endl;
                      //std::cout << "distfReduced" << std::endl << distfReduced << std::endl;
                    }

                    // Always use cell to the right of an interface to calculate upwinding for mu
                    if (imu < globalRgn.getUpper(4)-1)
                    {
                      idxUpwind[4] = idx[4] + 1;
                      distfIn.setPtr(distfInUpwindPtr, idxUpwind);
                      for (int nodeIndex = 0; nodeIndex < nodalStencilFixedVPar.size(); nodeIndex++)
                        distfUpwindReduced(nodeIndex) = distfInUpwindPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[nodeIndex]];
                      fUpperAvg = distfUpwindReduced.mean();
                    }

                    // Use cell to the right of the lower interface, i.e. the current cell
                    if (imu == globalRgn.getLower(4))
                      fLowerAvg = 0.0;
                    else
                      fLowerAvg = fOldAvg;


                    // Put in limiters on cell average?
                    if (fUpperAvg < 0.0)
                      fUpperAvg = 0.0;
                    if (fLowerAvg < 0.0)
                      fLowerAvg = 0.0;
                    
                    double alpha = 1.0/maxPosMu;
                    double fIncrement = (cellCentroid[4] + 0.5*grid.getDx(4))*fUpperAvg -
                      (cellCentroid[4] - 0.5*grid.getDx(4))*fLowerAvg;
                    double fNewAvg = fOldAvg + alpha*fIncrement;

                    fNewAvg = fOldAvg + fUpperAvg - fLowerAvg;
                    
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
                    for (int nodeIndex = 0; nodeIndex < nodalStencilFixedVPar.size(); nodeIndex++)
                    {
                      if (fOldAvg != 0.0)
                      {
                        distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[nodeIndex]] = 
                          (fNewAvg/fOldAvg)*distfInPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[nodeIndex]];
                      }
                      else
                      {
                        // Since f was initially 0 before drag step in this cell, set to a constant at these points
                        // to reach the desired density
                        //std::cout << "Setting new distribution function to a constant since initial density was zero." << std::endl;
                        //std::cout << "fNewAvg = " << fNewAvg << std::endl;
                        distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[nodeIndex]] = fNewAvg;
                      }
                    }
                    distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[0]] = 2*fNewAvg;
                    distfOutPtr[configNode + nodalStencilFixedMu[vSliceIndex] + nodalStencilFixedVPar[1]] = 0.0;*/
                    
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
  SOLPositivityDragSliceInMuUpdater::declareTypes()
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
  SOLPositivityDragSliceInMuUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLPositivityDragSliceInMuUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }

  bool
  SOLPositivityDragSliceInMuUpdater::sameSpecifiedCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList, const std::vector<int>& dimList)
  {
    for (int dimIndex = 0; dimIndex < dimList.size(); dimIndex++)
      if (std::fabs(nodeList(srcIndex,dimList[dimIndex])-nodeList(tarIndex,dimList[dimIndex])) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

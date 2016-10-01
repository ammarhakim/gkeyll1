/**
 * @file	LcSOLDesiredChargeDensity5DUpdater.cpp
 *
 * @brief	Adds source particles until a desired charge density has been reached
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLDesiredChargeDensity5DUpdater.h>

namespace Lucee
{
  const char *SOLDesiredChargeDensity5DUpdater::id = "SOLDesiredChargeDensity5D";

  SOLDesiredChargeDensity5DUpdater::SOLDesiredChargeDensity5DUpdater()
    : fnRef(-1)
  {
  }

  void
  SOLDesiredChargeDensity5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLDesiredChargeDensity5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLDesiredChargeDensity5DUpdater::readInput: Must specify element to use using 'basis3d'");
  }

  void
  SOLDesiredChargeDensity5DUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis5d->getNumNodes();
    std::vector<unsigned> zRef(nlocal), vRef(nlocal);

    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // Get a copy of the 5d nodal coordinates
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
  }

  Lucee::UpdaterStatus
  SOLDesiredChargeDensity5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // A 3D field containing the desired density at every node
    const Lucee::Field<5, double>& distfSource = this->getInp<Lucee::Field<5, double> >(0);
    // Density Delta is defined as selfDensity - otherDensity
    const Lucee::Field<3, double>& currDensityDelta = this->getInp<Lucee::Field<3, double> >(1);
    const Lucee::Field<3, double>& targetDensityDelta = this->getInp<Lucee::Field<3, double> >(2);
    // Distribution function to be scaled
    Lucee::Field<5, double>& distf = this->getOut<Lucee::Field<5, double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfSourcePtr = distfSource.createConstPtr();
    Lucee::ConstFieldPtr<double> currDensityDeltaPtr = currDensityDelta.createConstPtr();
    Lucee::ConstFieldPtr<double> targetDensityDeltaPtr = targetDensityDelta.createConstPtr();

    Lucee::FieldPtr<double> distfPtr = distf.createPtr(); // Output pointer

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int idx[5];
    double xc[5];

    // Loop over each node in configuration space
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[0] = ix;
          idx[1] = iy;
          idx[2] = iz;
          // Set density delta pointers to a position space cell
          currDensityDelta.setPtr(currDensityDeltaPtr, idx[0], idx[1], idx[2]);
          targetDensityDelta.setPtr(targetDensityDeltaPtr, idx[0], idx[1], idx[2]);

          // Loop over each 3d configuration space node
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // Figure out amount of distfSource to add
            double sourceDensity = targetDensityDeltaPtr[configNode] - currDensityDeltaPtr[configNode];
            // Only need to add more particles if sourceDensity is positive. If it is negative,
            // then more particles of the other species will be added in its own update
            if (sourceDensity > 0.0)
            {
              for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
              {
                idx[3] = iv;
                for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
                {
                  idx[4] = imu;
                  distfSource.setPtr(distfSourcePtr, idx);
                  distf.setPtr(distfPtr, idx);
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    distfPtr[nodalStencil[nodeIndex] + configNode] = distfPtr[nodalStencil[nodeIndex] + configNode] + 
                      sourceDensity*distfSourcePtr[nodalStencil[nodeIndex] + configNode];
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
  SOLDesiredChargeDensity5DUpdater::declareTypes()
  {
    // A 5D field containing a source distribution function to reach desired difference
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
    // Current value of density difference (e.g. between this species and another species)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // Target density difference (e.g. between this species and another species)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // Distribution function that will be modified to have the desired density given by
    // the input field
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLDesiredChargeDensity5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLDesiredChargeDensity5DUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

/**
 * @file	LcSOLDensityFloor5DUpdater.cpp
 *
 * @brief	Enforces a positive density floor in the domain by adding a source distribution function
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLDensityFloor5DUpdater.h>

namespace Lucee
{
  const char *SOLDensityFloor5DUpdater::id = "SOLDensityFloor5D";

  SOLDensityFloor5DUpdater::SOLDensityFloor5DUpdater()
    : fnRef(-1)
  {
  }

  void
  SOLDensityFloor5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLDensityFloor5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLDensityFloor5DUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLDensityFloor5DUpdater::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("densityFloor"))
      densityFloor = tbl.getNumber("densityFloor");
    else
      throw Lucee::Except("SOLDensityFloor5DUpdater::readInput: Must specify density floor using 'densityFloor'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLDensityFloor5DUpdater::initialize()
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

    // Compute matrix for v,mu integration
    int nlocal2d = nodalBasis2d->getNumNodes();
    Lucee::Matrix<double> tempMassMatrix(nlocal2d, nlocal2d);

    Eigen::MatrixXd integrationMatrix(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(tempMassMatrix);
    copyLuceeToEigen(tempMassMatrix, integrationMatrix);
    integrationMatrix *= grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));

    integrationVector = integrationMatrix.colwise().sum();
  }

  Lucee::UpdaterStatus
  SOLDensityFloor5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // A 3D field containing the desired density at every node
    const Lucee::Field<5, double>& distfSource = this->getInp<Lucee::Field<5, double> >(0);
    const Lucee::Field<3, double>& weightField = this->getInp<Lucee::Field<3, double> >(1);
    // Distribution function to be scaled
    Lucee::Field<5, double>& distf = this->getOut<Lucee::Field<5, double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfSourcePtr = distfSource.createConstPtr();
    Lucee::ConstFieldPtr<double> weightFieldPtr = weightField.createConstPtr();
    Lucee::FieldPtr<double> distfPtr = distf.createPtr(); // Output pointer

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int idx[5];
    double xc[5];
    Eigen::VectorXd distfReduced(nodalStencil.size());

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
          idx[3] = localRgn.getLower(3);
          idx[4] = localRgn.getLower(4);
          grid.setIndex(idx);
          // Set weighting factor for density calculation
          weightField.setPtr(weightFieldPtr, idx[0], idx[1], idx[2]);

          // Fill out absolute node list using local node coordinate list
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // Compute density of numerical solution at this configuration node
            double numericalDensity = 0.0;
            for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
            {
              idx[3] = iv;
              for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
              {
                idx[4] = imu;
                distf.setPtr(distfPtr, idx);
                for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  distfReduced(nodeIndex) = distfPtr[nodalStencil[nodeIndex] + configNode];
                numericalDensity += integrationVector.dot(distfReduced);
              }
            }

            // Scale numericalDensity so that its in the same units as the input targetDensityField
            numericalDensity = scaleFactor*weightFieldPtr[configNode]*numericalDensity;

            // If density is negative at this location, add enough of source to maintain density floor
            if (numericalDensity < densityFloor)
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
                      (densityFloor-numericalDensity)/densityFloor*distfSourcePtr[nodalStencil[nodeIndex] + configNode];
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
  SOLDensityFloor5DUpdater::declareTypes()
  {
    // A 5D field containing a source distribution function to offset negativity
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
    // A 3D field containing magnetic field weighting factoring for density calculation
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // Distribution function that will be modified to have the desired density given by
    // the input field
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLDensityFloor5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLDensityFloor5DUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

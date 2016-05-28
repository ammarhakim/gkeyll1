/**
 * @file	LcSOLMaxwellianAtNodeCalc.cpp
 *
 * @brief	Initializes a 5d maxwellian based on 3d number density, mean parallel velocity, and total temperature fields
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLMaxwellianAtNodeCalc.h>
#include <LcMathPhysConstants.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLMaxwellianAtNodeCalc::id = "SOLMaxwellianAtNodeCalc";

  SOLMaxwellianAtNodeCalc::SOLMaxwellianAtNodeCalc()
  {
  }

  void
  SOLMaxwellianAtNodeCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLMaxwellianAtNodeCalc::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLMaxwellianAtNodeCalc::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLMaxwellianAtNodeCalc::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("SOLMaxwellianAtNodeCalc::readInput: Must specify mass using 'speciesMass'");
  }

  void
  SOLMaxwellianAtNodeCalc::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis5d->getNumNodes();
    std::vector<unsigned> zRef(nlocal), vRef(nlocal);

    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // Get a copy of the 5d nodal coordinates
    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 5);
    nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
    nodeCoords = Eigen::MatrixXd(nlocal, 5);
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
  SOLMaxwellianAtNodeCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Input: number density
    const Lucee::Field<3, double>& numDensIn = this->getInp<Lucee::Field<3, double> >(0);
    // Input: <v> = n*u 
    const Lucee::Field<3, double>& mom1dir3In = this->getInp<Lucee::Field<3, double> >(1);
    // Input: temperature
    const Lucee::Field<3, double>& temperatureIn = this->getInp<Lucee::Field<3, double> >(2);
    // Input: magnetic field
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(3);
    // Output: distribution function
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);

    Lucee::ConstFieldPtr<double> numDensPtr = numDensIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1dir3Ptr = mom1dir3In.createConstPtr();
    Lucee::ConstFieldPtr<double> temperaturePtr = temperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();
    
    Lucee::FieldPtr<double> distfPtr = distfOut.createPtr();

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    double cellCentroid[5];
    int idx[5];
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<5> seq(localRgn);
    Lucee::Matrix<double> nodeCoordsLucee(nlocal5d, 5);

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      
      numDensIn.setPtr(numDensPtr, idx[0], idx[1], idx[2]);
      mom1dir3In.setPtr(mom1dir3Ptr, idx[0], idx[1], idx[2]);
      temperatureIn.setPtr(temperaturePtr, idx[0], idx[1], idx[2]);
      bFieldIn.setPtr(bFieldPtr, idx[0], idx[1], idx[2]);
      
      distfOut.setPtr(distfPtr, idx);
      
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);
      nodalBasis5d->setIndex(idx);

      // Get a copy of the 5d nodal coordinates
      nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
      copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

      // Loop over configuration space coordinates
      for (int configNode = 0; configNode < nlocal3d; configNode++)
      {
        double vThermSq = temperaturePtr[configNode]/speciesMass;
        double uVal = mom1dir3Ptr[configNode]/numDensPtr[configNode];
        // Loop over 5d nodes at this configuration space location
        for (int i = 0; i < nodalStencil.size(); i++)
        {
          double nodeIndex = configNode + nodalStencil[i];
          double vVal = nodeCoords(nodeIndex,3);
          double muVal = nodeCoords(nodeIndex,4);

          if (numDensPtr[configNode] == 0.0)
          {
            // This should not happen
            std::cout << "number density zero at point" << std::endl;
          }

          distfPtr[nodeIndex] = numDensPtr[configNode]/(2*PI*vThermSq*sqrt(2*PI*vThermSq))*exp(-(vVal-uVal)*(vVal-uVal)/(2*vThermSq))*
            exp(-muVal*bFieldPtr[configNode]/temperaturePtr[configNode]);
        }
      }
    }

    // Now scale so correct density is obtained
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
          bFieldIn.setPtr(bFieldPtr, idx[0], idx[1], idx[2]);
          numDensIn.setPtr(numDensPtr, idx[0], idx[1], idx[2]);

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
                distfOut.setPtr(distfPtr, idx);
                for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  distfReduced(nodeIndex) = distfPtr[nodalStencil[nodeIndex] + configNode];
                numericalDensity += integrationVector.dot(distfReduced);
              }
            }
            numericalDensity = numericalDensity*scaleFactor*bFieldPtr[configNode];

            // Loop through nodes a second time to scale to correct value
            for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
            {
              idx[3] = iv;
              for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
              {
                idx[4] = imu;
                distfOut.setPtr(distfPtr, idx);
                for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  distfPtr[nodalStencil[nodeIndex] + configNode] = (numDensPtr[configNode]/numericalDensity)*
                    distfPtr[nodalStencil[nodeIndex] + configNode];
              }
            }
          }
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLMaxwellianAtNodeCalc::declareTypes()
  {
    // Number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // <v> = u*n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Temperature (in joules)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Magnetic field profile
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: Maxwellian evaluated at nodes
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLMaxwellianAtNodeCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLMaxwellianAtNodeCalc::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

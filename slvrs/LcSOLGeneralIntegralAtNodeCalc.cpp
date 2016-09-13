/**
 * @file	LcSOLGeneralIntegralAtNodeCalc.cpp
 *
 * @brief	Compute an integral at each node in configuration space using a (v,mu) integration
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLGeneralIntegralAtNodeCalc.h>

namespace Lucee
{
  const char *SOLGeneralIntegralAtNodeCalc::id = "SOLGeneralIntegralAtNodeCalc";

  SOLGeneralIntegralAtNodeCalc::SOLGeneralIntegralAtNodeCalc()
  {
  }

  void
  SOLGeneralIntegralAtNodeCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLGeneralIntegralAtNodeCalc::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLGeneralIntegralAtNodeCalc::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLGeneralIntegralAtNodeCalc::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLGeneralIntegralAtNodeCalc::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis5d->getNumNodes();

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

    // Compute matrix for v,mu integration
    int nlocal2d = nodalBasis2d->getNumNodes();
    Lucee::Matrix<double> tempMassMatrix(nlocal2d, nlocal2d);
    integrationMatrix = Eigen::MatrixXd(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(tempMassMatrix);
    copyLuceeToEigen(tempMassMatrix, integrationMatrix);
    integrationMatrix *= grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));
  }

  Lucee::UpdaterStatus
  SOLGeneralIntegralAtNodeCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // Field to multiply distribution function by before integrating
    const Lucee::Field<5, double>& productIn = this->getInp<Lucee::Field<5, double> >(1);
    // Magnetic field in 3d
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(2);
    // Output field
    Lucee::Field<3, double>& integratedField = this->getOut<Lucee::Field<3, double> >(0);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> productInPtr = productIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();
    Lucee::FieldPtr<double> integratedFieldPtr = integratedField.createPtr(); // Output pointer

    integratedField = 0.0; // clear out current contents
    
    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int idx[5];

    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    Eigen::VectorXd distfReduced(nodalStencil.size());
    Eigen::VectorXd productReduced(nodalStencil.size());

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // Set input field pointers
      distfIn.setPtr(distfInPtr, idx);
      productIn.setPtr(productInPtr, idx);
      bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);
      // Set output field pointer
      integratedField.setPtr(integratedFieldPtr, idx[0], idx[1], idx[2]);

      // Loop over the  configuration space vertices in this cell
      for (int configNode = 0; configNode < nlocal3d; configNode++)
      {
        // At this particular configuration space vertix, copy all
        // nodes that occupy this location to a vector
        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
        {
          distfReduced(nodeIndex) = distfInPtr[nodalStencil[nodeIndex] + configNode];
          productReduced(nodeIndex) = productInPtr[nodalStencil[nodeIndex] + configNode];
        }

        // Accumulate results of integration
        integratedFieldPtr[configNode] += bFieldInPtr[configNode]*scaleFactor*
          distfReduced.dot(integrationMatrix*productReduced);
      }
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLGeneralIntegralAtNodeCalc::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: field to multiply with the distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: 3d field of magnetic field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: 3d field containing mean energy at each node
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  SOLGeneralIntegralAtNodeCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLGeneralIntegralAtNodeCalc::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}

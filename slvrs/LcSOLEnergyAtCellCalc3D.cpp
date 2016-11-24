/**
 * @file	LcSOLEnergyAtCellCalc3D.cpp
 *
 * @brief	Compute energy at each configuration space cell using a (x,y,z,v,mu) integration
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLEnergyAtCellCalc3D.h>

namespace Lucee
{
  const char *SOLEnergyAtCellCalc3D::id = "SOLEnergyAtCellCalc3D";

  SOLEnergyAtCellCalc3D::SOLEnergyAtCellCalc3D()
  {
  }

  void
  SOLEnergyAtCellCalc3D::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLEnergyAtCellCalc3D::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except("SOLEnergyAtCellCalc3D::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLEnergyAtCellCalc3D::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Get interpolation matrix for 3d integration
    // get volume quadrature data for 3d element
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    volWeights3d = std::vector<double>(nVolQuad3d);
    Lucee::Matrix<double> tempSurfQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempSurfCoords3d(nVolQuad3d, 3);

    nodalBasis3d->getGaussQuadData(tempSurfQuad3d, tempSurfCoords3d, volWeights3d);

    volQuad3d = Eigen::MatrixXd(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempSurfQuad3d, volQuad3d);
  }

  Lucee::UpdaterStatus
  SOLEnergyAtCellCalc3D::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Distribution function
    const Lucee::Field<3, double>& distfIn = this->getInp<Lucee::Field<3, double> >(0);
    // Hamiltonian
    const Lucee::Field<3, double>& hamilIn = this->getInp<Lucee::Field<3, double> >(1);
    // Total energy in a cell
    Lucee::Field<1, double>& energyField = this->getOut<Lucee::Field<1, double> >(0);
    
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilInPtr = hamilIn.createConstPtr();
    Lucee::FieldPtr<double> energyFieldPtr = energyField.createPtr(); // Output pointer
    
    energyField = 0.0; // clear out current contents

    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();

    int idx[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);

    Eigen::VectorXd distfVec(nlocal3d);
    Eigen::VectorXd hamilVec(nlocal3d);
    
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      
      energyField.setPtr(energyFieldPtr, idx[0]);
      distfIn.setPtr(distfInPtr, idx);
      hamilIn.setPtr(hamilInPtr, idx);

      // Fill out fields at nodes
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        distfVec(nodeIndex) = distfInPtr[nodeIndex];
        hamilVec(nodeIndex) = hamilInPtr[nodeIndex];
      }

      // Compute fields at quadrature points
      Eigen::VectorXd distfAtQuad = volQuad3d*distfVec;
      Eigen::VectorXd hamilAtQuad = volQuad3d*hamilVec;

      // Compute energy integral in cell
      double localTotalEnergy = 0.0;
      for (int quadIndex = 0; quadIndex < nVolQuad3d; quadIndex++)
        localTotalEnergy += volWeights3d[quadIndex]*distfAtQuad(quadIndex)*
          hamilAtQuad(quadIndex);

      // Convert integrated quantity to a cell-averaged physical quantity
      localTotalEnergy = scaleFactor*localTotalEnergy/grid.getDx(0);

      // Accumulate results of integration
      for (int nodeIndex = 0; nodeIndex < nlocal1d; nodeIndex++)
        energyFieldPtr[nodeIndex] += localTotalEnergy;
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLEnergyAtCellCalc3D::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: 1d field containing energy at each node
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  SOLEnergyAtCellCalc3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
}

/**
 * @file	LcSOLEnergyAtCellCalc.cpp
 *
 * @brief	Compute energy at each configuration space cell using a (x,y,z,v,mu) integration
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLEnergyAtCellCalc.h>

namespace Lucee
{
  const char *SOLEnergyAtCellCalc::id = "SOLEnergyAtCellCalc";

  SOLEnergyAtCellCalc::SOLEnergyAtCellCalc()
  {
  }

  void
  SOLEnergyAtCellCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLEnergyAtCellCalc::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLEnergyAtCellCalc::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLEnergyAtCellCalc::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Get interpolation matrix for 5d integration
    // get volume quadrature data for 5d element
    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    volWeights5d = std::vector<double>(nVolQuad5d);
    Lucee::Matrix<double> tempSurfQuad5d(nVolQuad5d, nlocal5d);
    Lucee::Matrix<double> tempSurfCoords5d(nVolQuad5d, 5);

    nodalBasis5d->getGaussQuadData(tempSurfQuad5d, tempSurfCoords5d, volWeights5d);

    volQuad5d = Eigen::MatrixXd(nVolQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfQuad5d, volQuad5d);
  }

  Lucee::UpdaterStatus
  SOLEnergyAtCellCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // bfield5d
    const Lucee::Field<5, double>& bFieldIn = this->getInp<Lucee::Field<5, double> >(1);
    // Hamiltonian
    const Lucee::Field<5, double>& hamilIn = this->getInp<Lucee::Field<5, double> >(2);
    // Output parallel velocity moments (0, 1) vs time
    Lucee::Field<3, double>& energyField = this->getOut<Lucee::Field<3, double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilInPtr = hamilIn.createConstPtr();
    Lucee::FieldPtr<double> energyFieldPtr = energyField.createPtr(); // Output pointer
    
    energyField = 0.0; // clear out current contents

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();

    int idx[5];
    Lucee::RowMajorSequencer<5> seq(localRgn);

    Eigen::VectorXd distfVec(nlocal5d);
    Eigen::VectorXd bFieldVec(nlocal5d);
    Eigen::VectorXd hamilVec(nlocal5d);
    
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      
      energyField.setPtr(energyFieldPtr, idx[0], idx[1], idx[2]);
      distfIn.setPtr(distfInPtr, idx);
      bFieldIn.setPtr(bFieldInPtr, idx);
      hamilIn.setPtr(hamilInPtr, idx);

      // Fill out fields at nodes
      for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
      {
        distfVec(nodeIndex) = distfInPtr[nodeIndex];
        hamilVec(nodeIndex) = hamilInPtr[nodeIndex];
        bFieldVec(nodeIndex) = bFieldInPtr[nodeIndex];
      }

      // Compute fields at quadrature points
      Eigen::VectorXd distfAtQuad = volQuad5d*distfVec;
      Eigen::VectorXd hamilAtQuad = volQuad5d*hamilVec;
      Eigen::VectorXd bFieldAtQuad = volQuad5d*bFieldVec;

      // Compute energy integral in cell
      double localTotalEnergy = 0.0;
      for (int quadIndex = 0; quadIndex < nVolQuad5d; quadIndex++)
        localTotalEnergy += volWeights5d[quadIndex]*distfAtQuad(quadIndex)*
          bFieldAtQuad(quadIndex)*hamilAtQuad(quadIndex);

      // Accumulate results of integration
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
        energyFieldPtr[nodeIndex] += scaleFactor*localTotalEnergy;
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLEnergyAtCellCalc::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: bfield5d
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Output: 3d field containing energy at each node
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  SOLEnergyAtCellCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
